"""Dry-run regression for an existing discrete run switching to TRANS."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SNAKEMAKE = os.environ.get("SNAKEMAKE") or shutil.which("snakemake")


@unittest.skipUnless(SNAKEMAKE, "Snakemake is not installed")
class WorkflowContractTests(unittest.TestCase):
    def test_missing_prepared_trans_model_reschedules_ref_pgs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            artifact = root / "reference/pgs_score_files/ldpred2/example"
            artifact.mkdir(parents=True)
            prefix = artifact / "ref-example"
            Path(str(prefix) + ".score.gz").touch()
            pcs = root / "reference/pc_score_files/TRANS/ref-TRANS-pcs.EUR.scale"
            pcs.parent.mkdir(parents=True)
            pcs.touch()
            prepared = root / "prepared.tsv"
            prepared.write_text(f"name\tmethod\tpath\tlabel\nexample\tldpred2\t{artifact}\tExample\n")
            genotype = root / "genotype.txt"
            genotype.touch()
            targets = root / "target.txt"
            targets.write_text(f"name path type indiv_report unrel\nsubject {genotype} 23andMe T NA\n")
            config = root / "config.json"
            config.write_text(json.dumps({
                "outdir": str(root), "resdir": str(root / "resources"), "refdir": "NA",
                "config_file": str(config), "target_list": str(targets), "gwas_list": "NA",
                "score_list": "NA", "prepared_score_list": str(prepared),
                "pgs_methods": [], "pgs_scaling": ["continuous"], "testing": "chr22",
            }))
            for suffix in ("-EUR.profiles", "-EUR.scale", "-TRANS.scale", "-TRANS.profiles"):
                Path(str(prefix) + suffix).touch()
            (root / "reference/pgs_score_files/ref_pgs.done").touch()
            env = dict(os.environ, CONDA_DEFAULT_ENV="genopred", GIT_DISCOVERY_ACROSS_FILESYSTEM="1")
            result = subprocess.run([
                SNAKEMAKE, "--snakefile", str(ROOT / "pipeline/Snakefile"),
                "--directory", str(root), "--configfile", str(ROOT / "pipeline/config.yaml"), str(config),
                "--rerun-triggers", "mtime", "--cores", "1", "--dry-run", "ref_pgs",
            ], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("Missing output files", result.stdout)
            self.assertIn("ref-example-TRANS.model.rds", result.stdout)

            # Once the reference outputs exist, the target profile itself is a
            # tracked output and schedules target scoring when missing.
            Path(str(prefix) + "-TRANS.model.rds").touch()
            checks = root / "reference/target_checks/subject"
            checks.mkdir(parents=True)
            (checks / "ancestry_reporter.done").touch()
            (checks / "pc_projection-TRANS.done").touch()
            target = root / "subject/pgs/TRANS/ldpred2/example/subject-example-TRANS.profiles"
            result = subprocess.run([
                SNAKEMAKE, "--snakefile", str(ROOT / "pipeline/Snakefile"),
                "--directory", str(root), "--configfile", str(ROOT / "pipeline/config.yaml"), str(config),
                "--rerun-triggers", "mtime", "--cores", "1", "--dry-run", str(target),
            ], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("rule target_pgs_i:", result.stdout)
            self.assertIn(str(target), result.stdout)


if __name__ == "__main__":
    unittest.main()
