import tempfile
import unittest
from pathlib import Path
from scrna_workflow.runner import digest, valid_checkpoint, DEPS, execute

class RunnerTests(unittest.TestCase):
    def test_artifact_tampering_invalidates_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a'; p.write_text('source')
            r={'status':'completed','artifact_hashes':{str(p):digest(p)}}
            self.assertTrue(valid_checkpoint(r))
            p.write_text('changed'); self.assertFalse(valid_checkpoint(r))
    def test_dependency_graph_acyclic(self):
        done=set()
        while len(done)<len(DEPS):
            ready={n for n,ds in DEPS.items() if n not in done and set(ds)<=done}
            self.assertTrue(ready); done.update(ready)
    def test_output_lock_prevents_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'.run.lock'; p.write_text('active')
            with self.assertRaisesRegex(RuntimeError,'locked'): execute({'output_dir':d})
            self.assertEqual(p.read_text(),'active')

if __name__=='__main__': unittest.main()
