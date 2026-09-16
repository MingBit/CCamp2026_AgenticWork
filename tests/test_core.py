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
