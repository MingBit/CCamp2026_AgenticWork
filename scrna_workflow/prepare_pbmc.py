"""Prepare public 10x PBMC3k raw counts with barcode-aligned tutorial annotations.

Downloads are deliberately separate; this command reads local public files only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import anndata as ad
import pandas as pd
import scipy.io
from scipy import sparse
import yaml

COUNT_URL='https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
REFERENCE_URL='https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad'
EXPECTED_COUNT_HASH='847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5'

def prepare(directory, output):
    d=Path(directory).resolve(); out=Path(output).resolve(); out.mkdir(parents=True,exist_ok=True)
    archive=d/'pbmc3k_filtered_gene_bc_matrices.tar.gz'
    reference=d/'pbmc3k_processed_reference.h5ad'
    sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    if sha(archive)!=EXPECTED_COUNT_HASH: raise ValueError('Public 10x archive hash differs from documented Scanpy tutorial hash.')
    with tarfile.open(archive,'r:gz') as t:
        def member(name):
            matches=[m for m in t.getmembers() if m.isfile() and m.name.endswith('/'+name)]
            if len(matches)!=1: raise ValueError('Expected exactly one '+name)
            return t.extractfile(matches[0])
        matrix=sparse.csr_matrix(scipy.io.mmread(member('matrix.mtx')).T)
        genes=pd.read_csv(member('genes.tsv'),sep='\t',header=None,names=['gene_ids','gene_symbols'])
        barcodes=pd.read_csv(member('barcodes.tsv'),sep='\t',header=None)[0].astype(str)
    a=ad.AnnData(matrix,obs=pd.DataFrame(index=pd.Index(barcodes,name='cell_id')),var=genes.set_index('gene_symbols'))
    a.var_names_make_unique()
    ref=ad.read_h5ad(reference)
    if 'louvain' not in ref.obs: raise ValueError('Reference does not have expected tutorial louvain annotations')
    if not ref.obs_names.is_unique or not ref.obs_names.isin(a.obs_names).all(): raise ValueError('Reference barcode mismatch')
    a.obs['reference_cell_type']=ref.obs['louvain'].astype(str).reindex(a.obs_names).fillna('unavailable')
    a.obs['sample']='pbmc3k'; a.obs['donor']='healthy_donor_1'; a.obs['condition']='healthy'
    a.uns['reference_annotation_provenance']='Scanpy processed PBMC3k tutorial labels, same biological data; not independently measured cell identity'
    a.uns['raw_count_provenance']=COUNT_URL
    path=out/'pbmc3k_counts_annotated.h5ad'; a.write_h5ad(path,compression='gzip')
    manifest={'count_url':COUNT_URL,'reference_url':REFERENCE_URL,'raw_archive_sha256':sha(archive),
              'reference_sha256':sha(reference),'prepared_sha256':sha(path),'cells':a.n_obs,'genes':a.n_vars,
              'reference_label_counts':a.obs['reference_cell_type'].value_counts().to_dict(),
              'reference_use':'benchmark only; not used to choose representation or clustering',
              'limitations':['One biological donor, one condition; condition-level DE unsupported.','Reference labels inferred from the same expression data, not independent truth.'],
              'marker_source':'https://satijalab.org/seurat/articles/pbmc3k_tutorial.html'}
    (out/'public_input_manifest.json').write_text(json.dumps(manifest,indent=2))
    config={'input_path':str(path),'output_dir':str(out/'run'),'matrix_kind':'counts','organism':'human','tissue':'peripheral blood mononuclear cells',
            'sample_column':'sample','donor_column':'donor','condition_column':'condition','seed':17,'workers':2,
            'n_neighbors':15,'n_pcs':30,'n_top_genes':2000,'resolutions':[0.4,0.8,1.2],
            'biological_question':'Recover broad PBMC populations and test exploratory cell-type marker expression; public software benchmark.',
            'markers':{'B cells':['MS4A1','CD79A','CD79B'], 'T cells':['CD3D','CD3E','IL7R','CCR7'],
                       'NK cells':['GNLY','NKG7','KLRD1','PRF1'], 'Monocytes':['LYZ','S100A8','S100A9','FCGR3A','MS4A7','LST1'],
                       'Dendritic cells':['FCER1A','CST3','CD1C'], 'Platelets':['PPBP','PF4']}}
    (out/'config.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
    print(json.dumps(manifest,indent=2))
    return out/'config.yaml'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--downloads',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();prepare(args.downloads,args.output)
