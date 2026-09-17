"""Dependency scheduler and content-addressed checkpoints. No external data transfer."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import traceback
import yaml

DEPS = {'inspection': [], 'qc': ['inspection'], 'representation': ['qc'],
        'graph': ['representation'], 'clustering': ['graph'],
        'regulon': ['clustering'], 'discovery': ['clustering'],
        'validation': ['qc', 'graph', 'clustering', 'regulon', 'discovery'],
        'report': ['validation']}
REVIEW_TASKS = {'validation', 'report'}
# Input files hashed into the run fingerprint; path keys are also resolved relative to the config file.
INPUT_KEYS = ('input_path','metadata_path','tf_list','motif_evidence_path',
              'scenicplus_fragments_path','scenicplus_ctx_db_path','scenicplus_dem_db_path',
              'scenicplus_motif_annotations_path','scenicplus_blacklist_path','scenicplus_tss_annotation_path',
              'scenicplus_genome_annotation_path','scenicplus_chromsizes_path','scenicplus_mallet_path')
PATH_KEYS = INPUT_KEYS + ('output_dir','scenicplus_python','scenicplus_work_dir','scenicplus_temp_dir')
LARGE_INPUT_BYTES = 256 * 1024 * 1024

def digest(path):
    p = Path(path)
    if p.is_dir():
        return hashlib.sha256(json.dumps([(str(f.relative_to(p)), digest(f)) for f in sorted(p.rglob('*')) if f.is_file()]).encode()).hexdigest()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def input_digest(path):
    """SHA-256 of an input; large immutable files are cached by resolved path, size and mtime."""
    p = Path(path).resolve()
    if not p.is_file() or p.stat().st_size < LARGE_INPUT_BYTES:
        return digest(p)
    cache = Path(os.environ.get('SCRNA_WORKFLOW_DIGEST_CACHE', Path.home() / '.cache' / 'scrna_workflow' / 'input_digests.json'))
    stat = p.stat(); key = f'{p}|{stat.st_size}|{stat.st_mtime_ns}'
    try: known = json.loads(cache.read_text())
    except (OSError, ValueError): known = {}
    if key not in known:
        known[key] = digest(p)
        try:
            cache.parent.mkdir(parents=True, exist_ok=True); write_json(cache, known)
        except OSError: pass
    return known[key]

def write_json(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, default=str) + '\n')
    tmp.replace(path)

def valid_checkpoint(result):
    hashes = result.get('artifact_hashes', {})
    return result.get('status') in ('success', 'completed', 'passed', 'skipped', 'inconclusive') and all(Path(p).is_file() and digest(p) == h for p,h in hashes.items())

def execute(config, resume=False):
    root = Path(config.get('output_dir', 'results')).expanduser().resolve()
    source = Path(config['input_path']).expanduser().resolve() if config.get('input_path') else None
    if source and source.is_dir() and (source == root or source in root.parents):
        raise ValueError('Output directory must be outside the input directory to preserve immutable input hashing.')
    root.mkdir(parents=True, exist_ok=True)
    cache = root / '.cache'
    cache.mkdir(exist_ok=True)
    os.environ.setdefault('NUMBA_CACHE_DIR', str(cache / 'numba'))
    os.environ.setdefault('MPLCONFIGDIR', str(cache / 'matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(cache))
    from . import core, downstream
    # Single orchestrator owns run state; each worker exclusively owns its task directory.
    lock = root / '.run.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(f'Output workspace locked: {lock}. Confirm no process is active before removing a stale lock.')
    os.write(fd, str(os.getpid()).encode()); os.close(fd)
    try:
        return _execute(config, root, resume, core, downstream)
    finally:
        lock.unlink(missing_ok=True)

def _execute(config, root, resume, core, downstream):
    versions = {}
    for name in ('numpy','scipy','pandas','scanpy','anndata','scikit-learn','matplotlib','igraph','leidenalg','PyYAML','pydeseq2','numba','umap-learn','statsmodels','networkx','h5py'):
        try: versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: versions[name] = None
    inputs = {}
    for key in INPUT_KEYS:
        value = config.get(key)
        if value and Path(value).exists():
            inputs[key] = {'path': str(Path(value).resolve()), 'sha256': input_digest(value)}
    source_hash = hashlib.sha256(json.dumps([(p.name,digest(p)) for p in sorted(Path(__file__).parent.glob('*.py'))]).encode()).hexdigest()
    fingerprint = hashlib.sha256(json.dumps({'config':config,'inputs':inputs,'source':source_hash,'versions':versions,'python':sys.version,'platform':platform.platform()},sort_keys=True).encode()).hexdigest()
    state_path = root / 'run_state.json'
    old = json.loads(state_path.read_text()) if state_path.exists() and resume else {}
    if old and old.get('fingerprint') != fingerprint:
        raise ValueError('Resume refused: inputs, code, parameters, or environment changed. Use a new output directory.')
    if state_path.exists() and not resume:
        raise ValueError('Output already contains a run; use --resume or a new output directory.')
    resources = {'logical_cpus': os.cpu_count(), 'gpu_required': False, 'workers': config.get('workers',2)}
    try: resources['physical_memory_bytes'] = os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE')
    except (ValueError, OSError): resources['physical_memory_bytes'] = None
    state = {'schema_version':1, 'fingerprint':fingerprint,'config':config,'inputs':inputs,
             'versions':versions,'resources':resources,'python':sys.version,'platform':platform.platform(),
             'source_sha256':source_hash,'seed':config.get('seed',0),'dependencies':DEPS,
             'tasks':old.get('tasks',{}),'started_at':old.get('started_at',datetime.now(timezone.utc).isoformat()),
             'decisions':['No automatic biological annotation without supplied marker evidence.',
                          'No repeated scientific failure: automatic retries limited to transient OSError.',
                          'Validation requests are recorded; revisions require a specific changed assumption and at most two cycles per issue.',
                          'Input data are never sent to external services.']}
    write_json(root/'manifest.json', {k:v for k,v in state.items() if k!='tasks'})
    functions = {'inspection':core.inspect_data,'qc':core.qc,'representation':core.representation,
                 'graph':core.graph,'clustering':core.clustering,**{n:getattr(downstream,n) for n in ('regulon','discovery','validation','report')}}
    artifacts = {}
    pending = set(DEPS)
    changed = set()
    log_path = root/'events.jsonl'
    def event(data):
        with log_path.open('a') as f: f.write(json.dumps({'time':datetime.now(timezone.utc).isoformat(),**data},default=str)+'\n')
    def call(name):
        task_dir = root/name; task_dir.mkdir(exist_ok=True)
        ctx = {'config':config,'root':root,'artifacts':dict(artifacts),'seed':config.get('seed',0)}
        for attempt in range(min(max(int(config.get('transient_retries',1)),0),2)+1):
            try:
                result = functions[name](ctx)
                required = {'status','input_references','outputs','metrics','warnings','recommended_next_actions'}
                if not required.issubset(result): raise ValueError(f'Missing result fields: {required-set(result)}')
                result['attempts'] = attempt+1
                result['artifact_hashes'] = {str(p):digest(p) for p in result['outputs'].values() if isinstance(p,str) and Path(p).is_file()}
                write_json(task_dir/'result.json',result)
                return result
            except Exception as exc:
                if isinstance(exc,OSError) and attempt < min(max(int(config.get('transient_retries',1)),0),2): continue
                result = {'status':'failed','input_references':DEPS[name],'outputs':{},'metrics':{},
                          'warnings':[f'{type(exc).__name__}: {exc}'],'recommended_next_actions':['Inspect traceback and correct cause before a deliberate rerun.'],
                          'traceback':traceback.format_exc(),'attempts':attempt+1}
                write_json(task_dir/'result.json',result)
                return result
    while pending:
        ready = sorted(n for n in pending if all(d in artifacts for d in DEPS[n]))
        if not ready: raise RuntimeError('Dependency cycle')
        runnable = []
        for name in ready:
            previous = state['tasks'].get(name)
            if previous and not any(d in changed for d in DEPS[name]) and valid_checkpoint(previous):
                artifacts[name] = previous; pending.remove(name); event({'task':name,'event':'checkpoint_reused'}); continue
            if previous and previous.get('status') == 'failed':
                artifacts[name]=previous; pending.remove(name); event({'task':name,'event':'failure_preserved','reason':'No changed cause supplied; start a new run after correction.'}); continue
            blocked = [d for d in DEPS[name] if artifacts[d]['status'] in ('failed','blocked')]
            if blocked and name not in REVIEW_TASKS:
                artifacts[name] = {'status':'blocked','input_references':blocked,'outputs':{},'metrics':{},'warnings':['Required input unavailable'],'recommended_next_actions':['Resolve upstream blockers.']}
                state['tasks'][name] = artifacts[name]; pending.remove(name)
            else: changed.add(name); runnable.append(name); event({'task':name,'event':'call','function':functions[name].__name__,'dependencies':DEPS[name]})
        with ThreadPoolExecutor(max_workers=max(1,min(int(config.get('workers',2)),4))) as pool:
            futures = {pool.submit(call,n):n for n in runnable}
            for future in as_completed(futures):
                name=futures[future]; result=future.result(); artifacts[name]=result
                state['tasks'][name]=result; pending.remove(name)
                event({'task':name,'event':'returned','status':result['status'],'warnings':result['warnings']})
                write_json(state_path,state)
        write_json(state_path,state)
    state['completed_at']=datetime.now(timezone.utc).isoformat()
    state['status']='blocked' if any(r['status'] in ('blocked','failed') for r in artifacts.values()) else 'completed'
    write_json(state_path,state)
    print(json.dumps({'status':state['status'],'output_dir':str(root),'tasks':{k:v['status'] for k,v in artifacts.items()}},indent=2))
    return 2 if state['status']=='blocked' else 0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default=None)
    p.add_argument('--input',dest='input_path'); p.add_argument('--metadata',dest='metadata_path')
    p.add_argument('--output',dest='output_dir'); p.add_argument('--resume',action='store_true')
    p.add_argument('--seed',type=int); p.add_argument('--workers',type=int)
    args=p.parse_args(); config={}
    if args.config:
        config=yaml.safe_load(Path(args.config).read_text()) or {}
        base=Path(args.config).resolve().parent
        for key in PATH_KEYS:
            if config.get(key): config[key]=str((base/Path(config[key]).expanduser()).resolve())
    for k in ('input_path','metadata_path','output_dir','seed','workers'):
        if getattr(args,k) is not None: config[k]=getattr(args,k)
    for key in PATH_KEYS:
        if config.get(key): config[key]=str(Path(config[key]).expanduser().resolve())
    try: sys.exit(execute(config,args.resume))
    except Exception as exc: print(f'ERROR: {exc}',file=sys.stderr); sys.exit(2)
