"""Small end-to-end regression test, including count and checkpoint integrity."""
import json
from pathlib import Path
import tempfile
import unittest
import yaml
from scrna_workflow.runner import execute, digest
from scrna_workflow.synthetic import generate

class IntegrationTests(unittest.TestCase):
    def test_synthetic_run_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            config=yaml.safe_load(generate(d).read_text())
            self.assertEqual(execute(config),0)
            root=Path(config['output_dir'])
            state=json.loads((root/'run_state.json').read_text())
            for task in ('inspection','qc','representation','graph','clustering','regulon','discovery','validation','report'):
                self.assertIn(state['tasks'][task]['status'],('completed','success','passed'))
            import anndata as ad
            from scipy import sparse
            original=ad.read_h5ad(config['input_path'])
            final=ad.read_h5ad(state['tasks']['clustering']['outputs']['dataset'])
            aligned=original[final.obs_names,final.var_names]
            self.assertEqual((sparse.csr_matrix(aligned.X)-sparse.csr_matrix(final.layers['counts'])).nnz,0)
            self.assertTrue(final.obs_names.is_unique)
            self.assertTrue(final.var_names.is_unique)
            for matrix in final.obsp.values(): self.assertEqual(matrix.shape,(final.n_obs,final.n_obs))
            before={p:digest(p) for r in state['tasks'].values() for p in r.get('artifact_hashes',{})}
            self.assertEqual(execute(config,resume=True),0)
            self.assertEqual(before,{p:digest(p) for p in before})
            self.assertIn('checkpoint_reused',(root/'events.jsonl').read_text())
            config['seed']=999
            with self.assertRaisesRegex(ValueError,'Resume refused'): execute(config,resume=True)

if __name__=='__main__': unittest.main()
