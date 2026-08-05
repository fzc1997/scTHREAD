#!/usr/bin/env python3
"""Unit tests for shipped NAR figure data loaders (real registry + figdata paths).

Run: python3 test_nar_loaders.py
"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_nar_bio as R  # noqa: E402
import nar_style as S  # noqa: E402


class TestCatalogHeadline(unittest.TestCase):
    def test_catalog_headline_positive_and_live_shape(self):
        h = R.catalog_headline()
        self.assertIsInstance(h["n_runs"], int)
        self.assertIsInstance(h["n_studies"], int)
        self.assertIsInstance(h["n_cells"], int)
        self.assertGreater(h["n_runs"], 0)
        self.assertGreater(h["n_studies"], 0)
        self.assertGreater(h["n_cells"], 100_000)
        # live baseline at goal time (must track registry, not hard-code forever if registry grows —
        # but assert consistency: runs >= studies, cells >> runs)
        self.assertGreaterEqual(h["n_runs"], h["n_studies"])
        self.assertGreater(h["n_cells"], h["n_runs"] * 100)
        # human+mouse partition
        self.assertGreaterEqual(h["n_human"] + h["n_mouse"], h["n_runs"] - 20)

    def test_catalog_matches_pinned_scope_recompute(self):
        """Drive the real entry point and recompute independently from the same source."""
        h = R.catalog_headline()
        done = R.registry_done()
        done["isoquant_cells"] = pd.to_numeric(done.get("isoquant_cells"), errors="coerce")
        self.assertEqual(h["n_runs"], int(len(done)))
        self.assertEqual(h["n_studies"], int(done["gse"].nunique()))
        self.assertEqual(h["n_cells"], int(done["isoquant_cells"].fillna(0).sum()))

    def test_default_scope_is_the_manuscript_snapshot(self):
        """Fig.1 must not drift back onto the growing live registry.

        The manuscript denominators are frozen; rendering from sample_registry.tsv
        silently produced 478/33/941,733 against a text that said 434/30/850,938.
        """
        self.assertEqual(R.CATALOG_SCOPE, "snapshot")
        h = R.catalog_headline()
        self.assertEqual(
            (h["n_runs"], h["n_studies"], h["n_cells"]), (434, 30, 850_938)
        )
        self.assertEqual((h["n_human"], h["n_mouse"]), (110, 324))
        self.assertEqual((h["n_ont"], h["n_pacbio"]), (384, 50))


class TestLoadCatalogBio(unittest.TestCase):
    def test_load_catalog_bio_all_done_gse(self):
        df = R.load_catalog_bio()
        self.assertGreater(len(df), 10)
        self.assertIn("gse", df.columns)
        self.assertIn("n_cells", df.columns)
        self.assertIn("system", df.columns)
        self.assertTrue((df["n_cells"] >= 0).all())
        h = R.catalog_headline()
        self.assertEqual(int(df["n_cells"].sum()), h["n_cells"])
        self.assertEqual(int(df["gse"].nunique()), h["n_studies"])


class TestFigdataIntegrity(unittest.TestCase):
    def test_cohort_headline_labeled_subset(self):
        ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())
        self.assertEqual(ch["n_studies"], 15)
        self.assertEqual(ch["n_runs"], 67)
        self.assertEqual(ch["n_annotated_cells"], 452_213)

    def test_marrow_reproducibility_rows(self):
        m = pd.read_csv(S.FIGDATA / "fig4b_marrow_reproducibility.tsv", sep="\t")
        self.assertGreaterEqual(len(m), 4)
        self.assertIn("novel_frac_myeloma", m.columns)
        self.assertIn("novel_frac_ccus", m.columns)

    def test_three_axis_tables_exist(self):
        for name, path in [
            ("ASE", S.F2DATA / "ase_interaction.tsv"),
            ("DIU", S.F2DATA / "diu_celltype.tsv"),
            ("APA", S.F2DATA / "apa_celltype.tsv"),
        ]:
            self.assertTrue(path.exists(), f"missing {path}")
            df = pd.read_csv(path, sep="\t")
            self.assertGreater(len(df), 1000, name)
            self.assertIn("sig", df.columns)


if __name__ == "__main__":
    unittest.main(verbosity=2)
