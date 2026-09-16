"""Make synthetic counts for software testing only; no biological interpretation."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import yaml

def generate(output):
    out=Path(output).resolve(); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(17)
    n,g=240,120
    groups=np.arange(n)%3
    rates=np.full((n,g),0.35)
    for k in range(3): rates[groups==k,k*20:(k+1)*20]=4
    counts=rng.poisson(rates*rng.lognormal(0,0.2,(n,1))).astype(np.int32)
    obs=pd.DataFrame({'sample':[f's{i//40}' for i in range(n)],'donor':[f'd{i//40}' for i in range(n)],
                      'condition':['A' if i//40<3 else 'B' for i in range(n)],'synthetic_group':groups.astype(str)},index=[f'cell{i:04}' for i in range(n)])
    genes=[f'G{i:03}' for i in range(g)]; genes[-3:]=['MT-SYN1','MT-SYN2','MT-SYN3']
    a=ad.AnnData(sparse.csr_matrix(counts),obs=obs,var=pd.DataFrame(index=genes))
    a.uns['synthetic_software_test']=True
    a.write_h5ad(out/'synthetic.h5ad')
    (out/'tf_list.txt').write_text('G000\nG020\nG040\n')
    config={'input_path':str(out/'synthetic.h5ad'),'output_dir':str(out/'run'),'seed':17,'workers':2,'matrix_kind':'counts',
            'sample_column':'sample','donor_column':'donor','condition_column':'condition','primary_comparison':['A','B'],'organism':'synthetic',
            'mitochondrial_prefix':'MT-','n_neighbors':12,'n_pcs':15,'n_top_genes':100,'resolutions':[0.4,0.8],
            'tf_list':str(out/'tf_list.txt'),'tf_resource_organism':'synthetic',
            'synthetic':True,'biological_question':'Software integration test only','markers':{}}
    (out/'config.yaml').write_text(yaml.safe_dump(config))
    return out/'config.yaml'

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True); print(generate(p.parse_args().output))
