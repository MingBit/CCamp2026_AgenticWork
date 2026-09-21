import os
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/scrna-test-numba')
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scrna_workflow.core import inspect_data, _countlike
from scrna_workflow.modalities import detect_modalities
from scrna_workflow.tools.gene_filter import remove_ribosomal_genes
from scrna_workflow.tools.annotation import annotate_clusters


class CoreTests(unittest.TestCase):
    def fixture(self, root):
        a = ad.AnnData(sparse.csr_matrix([[1,0,2],[0,3,1],[2,2,1]]),
                       obs=pd.DataFrame(index=['a','b','c']), var=pd.DataFrame(index=['g1','g2','g3']))
        p = Path(root) / 'input.h5ad'; a.write_h5ad(p)
        return a, {'root':Path(root)/'output','config':{'input_path':str(p)},'artifacts':{},'seed':1}

    def test_noncounts_not_misrepresented(self):
        self.assertFalse(_countlike(sparse.csr_matrix([[1.2,0],[0,2]])))
        self.assertFalse(_countlike(np.array([[1,-1],[0,2]])))

    def test_source_and_counts_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            source, ctx = self.fixture(d)
            r = inspect_data(ctx)
            observed = ad.read_h5ad(r['outputs']['dataset'])
            self.assertEqual((observed.layers['counts'] != source.X).nnz, 0)
            self.assertTrue((Path(d)/'output/inspection/source/input.h5ad').exists())

    def test_metadata_alignment_required(self):
        with tempfile.TemporaryDirectory() as d:
            _, ctx = self.fixture(d)
            p=Path(d)/'meta.csv'; pd.DataFrame({'donor':['x','y']},index=['a','b']).to_csv(p)
            ctx['config']['metadata_path']=str(p)
            with self.assertRaisesRegex(ValueError,'coverage'):
                inspect_data(ctx)

    def test_declared_counts_must_be_integer(self):
        with tempfile.TemporaryDirectory() as d:
            a, ctx = self.fixture(d)
            a.X = a.X.astype(float) * .3; a.write_h5ad(ctx['config']['input_path'])
            ctx['config']['matrix_kind']='counts'
            with self.assertRaisesRegex(ValueError,'noninteger'):
                inspect_data(ctx)

    def test_rpl_rps_filter_preserves_other_count_columns(self):
        a = ad.AnnData(
            sparse.csr_matrix([[1, 2, 3, 4, 5], [5, 6, 7, 8, 9]]),
            var=pd.DataFrame(index=['RPLP0', 'RPS3', 'CD3D', 'MS4A1', 'LST1']),
        )
        a.layers['counts'] = a.X.copy()
        filtered, removed = remove_ribosomal_genes(a)
        self.assertEqual(filtered.var_names.tolist(), ['CD3D', 'MS4A1', 'LST1'])
        self.assertEqual(removed['feature_id'].tolist(), ['RPLP0', 'RPS3'])
        self.assertEqual((filtered.layers['counts'] != a[:, ['CD3D', 'MS4A1', 'LST1']].layers['counts']).nnz, 0)
        self.assertEqual(a.n_vars, 5)

    def test_marker_annotation_requires_multiple_supporting_genes(self):
        a = ad.AnnData(np.ones((4, 4)),
                       obs=pd.DataFrame({'cluster': pd.Categorical(['0','0','1','1'])}),
                       var=pd.DataFrame(index=['CD3D','CD3E','MS4A1','CD79A']))
        markers = pd.DataFrame({'group':['0','0','1','1'],
                                'names':['CD3D','CD3E','MS4A1','CD79A'],
                                'logfoldchanges':[2.,2.,2.,2.]})
        result = annotate_clusters(a, markers, {'T cells':['CD3D','CD3E'],
                                                'B cells':['MS4A1','CD79A']})
        self.assertEqual([row['label'] for row in result], ['T cells','B cells'])
        self.assertTrue((a.obs['annotation_confidence'] > 0).all())

    def test_modality_detection_defaults_h5ad_to_rna(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'input.h5ad'
            ad.AnnData(sparse.csr_matrix([[1, 0], [0, 1]])).write_h5ad(path)
            self.assertEqual(detect_modalities(path), ['rna'])

    def test_modality_detection_reads_h5ad_feature_types(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'input.h5ad'
            data = ad.AnnData(sparse.csr_matrix([[1, 0], [0, 1]]),
                              var=pd.DataFrame({'feature_types': ['Peaks', 'Antibody Capture']}, index=['p1', 'a1']))
            data.write_h5ad(path)
            self.assertEqual(detect_modalities(path), ['atac', 'protein'])
