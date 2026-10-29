"""Focused checks without importing the Snakemake workflow."""
import ast
import os
from pathlib import Path
import re
import tempfile
import unittest

source = Path(__file__).resolve().parents[1] / "rules" / "dependencies.smk"
text = source.read_text().split('def check_target_paths(', 1)[1].split('\n# Create function', 1)[0]
tree = ast.parse('def check_target_paths(' + text)
function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "check_target_paths")
namespace = {"os": os}
exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
check_target_paths = namespace["check_target_paths"]
identifier_source = source.read_text().split('def check_identifiers(', 1)[1].split('\n###', 1)[0]
identifiers = ast.parse('def check_identifiers(' + identifier_source)
identifier_function = next(node for node in identifiers.body if isinstance(node, ast.FunctionDef))
class PandasStub:
    @staticmethod
    def isna(value):
        return value is None
identifier_namespace = {"re": re, "pd": PandasStub}
exec(compile(ast.Module(body=[identifier_function], type_ignores=[]), str(source), "exec"), identifier_namespace)
check_identifiers = identifier_namespace["check_identifiers"]


class Rows:
    def __init__(self, rows):
        self.rows = rows

    def iterrows(self):
        return enumerate(self.rows)


class Names:
    columns = ['name']

    def __init__(self, values):
        self.values = values

    def __getitem__(self, column):
        return self.values


class TargetPathTests(unittest.TestCase):
    def test_identifiers_cannot_interpolate_shell_or_traverse(self):
        check_identifiers(Names(['sample_1', 'GCST123.1']), 'name', 'test')
        for invalid in ('../escape', '..', '$(touch bad)', 'sample;bad', ''):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                check_identifiers(Names([invalid]), 'name', 'test')

    def test_all_rows_and_companions_are_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            for extension in (".bed", ".bim", ".fam"):
                Path(f"{path}.chr1{extension}").touch()
            valid = {"type": "plink1", "path": str(path)}
            missing = {"type": "23andMe", "path": str(path / "missing")}
            self.assertEqual(check_target_paths(Rows([valid]), "1"), [])
            with self.assertRaises(FileNotFoundError):
                check_target_paths(Rows([valid, missing]), "1")
            Path(f"{path}.chr1.bim").unlink()
            with self.assertRaises(FileNotFoundError):
                check_target_paths(Rows([valid]), "1")

    def test_bgen_requires_index_and_all_chromosomes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            Path(f"{path}.sample").touch()
            for chrom in (1, 2):
                Path(f"{path}.chr{chrom}.bgen").touch()
                Path(f"{path}.chr{chrom}.bgen.bgi").touch()
            row = {"type": "bgen", "path": str(path)}
            self.assertEqual(check_target_paths(Rows([row]), range(1, 3)), [])
            Path(f"{path}.chr2.bgen.bgi").unlink()
            with self.assertRaises(FileNotFoundError):
                check_target_paths(Rows([row]), range(1, 3))


if __name__ == "__main__":
    unittest.main()
