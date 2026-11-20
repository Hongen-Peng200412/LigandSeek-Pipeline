"""验证同一规范身份能穿过旧 Match 得分轴、人工改选和独立 Build 入口."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from ligandseek.identity import canonical_smiles, normalize_identities
from ligandseek.pipeline import run_sample


class IdentityTests(unittest.TestCase):
    """核对化学身份归并及历史得分文件的索引关系."""

    def test_canonical_smiles(self):
        """显式立体、同位素和显式氢经过多次规范化后仍得到同一身份."""
        for left, right in (("C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O"), ("F/C=C/F", r"F/C=C\F"),
                            ("[13CH3]CO", "CCO"), ("[2H]O[2H]", "O")):
            self.assertEqual(canonical_smiles(left), canonical_smiles(right))
            self.assertEqual(canonical_smiles(canonical_smiles(left)), canonical_smiles(left))

    def test_stored_axis_and_manual_build(self):
        """复用乱序且含重复身份的旧 Match 文件时, 规范化的手选身份独立交给 Build."""
        left, right, model = "C[C@H](F)Cl", "C[C@@H](F)Cl", "CC(F)Cl"
        with TemporaryDirectory() as scratch:
            root = Path(scratch)
            sample = {"sample_id": "demo", "experimental_map": "exp.mrc", "simulated_map": "sim.mrc", "receptor_cif": "rec.cif",
                      "identities": [{"smiles": text, "type_tag": "small_molecule"} for text in ("CCO", left, right)]}
            stored = {"identity_smiles": [right, "OCC", left], "candidates": [
                {"blob_id": 0, "center_xyz_A": [1, 2, 3], "selected": True, "probability_mean": 0.9,
                 "foreground_probability": 0.7, "identity_scores": [0.5, 0.7, 0.2]}]}
            for name, data in (("sample", sample), ("match", stored), ("decisions", [{"blob_id": 0, "action": "dock", "smiles": right}])):
                (root / (name + ".json")).write_text(json.dumps(data))
            with patch("ligandseek.pipeline.execute", return_value=[]) as execute:
                run_sample(root / "sample.json", root / "out", models=root, command="build",
                           candidates=root / "match.json", decisions=root / "decisions.json")
            plan = json.loads((root / "out/plan.json").read_text())
            self.assertEqual(plan["identity_smiles"], ["CCO", model])
            region, = plan["candidates"]
            self.assertEqual(region["identity_scores"], [0.7, 0.5])
            self.assertEqual(region["predicted_smiles"], "CCO")
            self.assertEqual(region["manual_input_smiles"], right)
            self.assertEqual(region["chosen_smiles"], model)
            self.assertEqual(execute.call_args.args[2][0]["chosen_smiles"], model)

    def test_alias_keeps_supplied_type_and_ccd_id(self):
        """同一身份的后续写法可以补足公开的化学类别、组分和向量, 已有字段优先."""
        sample = {"identities": [
            {"smiles": "OCC"},
            {"smiles": "CCO", "type_tag": "small_molecule", "ccd_id": "EOH", "embedding": "e.npz"},
            {"smiles": "OCC", "embedding": "late.npz"},
        ]}
        normalize_identities(sample, [], None, None)
        self.assertEqual(sample["identities"], [{"smiles": "CCO", "type_tag": "small_molecule",
                                                  "ccd_id": "EOH", "embedding": "e.npz"}])
