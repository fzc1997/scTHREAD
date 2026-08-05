#!/usr/bin/env python3
"""Promote the frozen 9,999-permutation screens into manuscript-facing tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
P0 = ROOT / "tables/p0_biological_unit_rerun"
SUMMARY = ROOT / "tables/Table_S_three_axis_summary.tsv"
PTPRC = ROOT / "tables/Table_S_PTPRC_multievidence.tsv"
INVENTORY = ROOT / "tables/NAR_FEATURED_INVENTORY.tsv"
GID = "ENSG00000081237"


def load(axis: str) -> pd.DataFrame:
    path = P0 / f"{axis}_observed_9999.tsv"
    frame = pd.read_csv(path, sep="\t")
    if frame.gene.duplicated().any() or not path.with_suffix(
        path.suffix + ".manifest.json"
    ).is_file():
        raise RuntimeError(f"invalid 9,999-permutation source: {path}")
    return frame


def one(frame: pd.DataFrame, gene: str, axis: str) -> pd.Series:
    rows = frame[frame.gene.astype(str).eq(gene)]
    if len(rows) != 1:
        raise RuntimeError(f"expected one {axis} row for {gene}; found {len(rows)}")
    return rows.iloc[0]


def main() -> None:
    frames = {axis: load(axis) for axis in ("ase", "diu", "apa")}
    summary_rows = []
    for axis in ("ase", "diu", "apa"):
        frame = frames[axis]
        summary_rows.append(
            {
                "axis": axis.upper(),
                "n_genes_tested": len(frame),
                "n_sig_fdr05_effect_gate": int(frame.sig.astype(bool).sum()),
                "frac_raw_p_lt_0.05": float(frame.pval.lt(0.05).mean()),
                "source": str(P0 / f"{axis}_observed_9999.tsv"),
            }
        )
    pd.DataFrame(summary_rows).to_csv(SUMMARY, sep="\t", index=False)

    diu = one(frames["diu"], GID, "DIU")
    apa = one(frames["apa"], GID, "APA")
    ase = one(frames["ase"], "PTPRC", "ASE")
    pd.DataFrame(
        [
            {
                "gene_symbol": "PTPRC",
                "ensembl_id": GID,
                "diu_sig": bool(diu.sig),
                "diu_qval": float(diu.qval),
                "diu_effect": float(diu.effect_equal_donor),
                "diu_n_iso": int(diu.n_features),
                "apa_sig": bool(apa.sig),
                "apa_qval": float(apa.qval),
                "apa_effect": float(apa.effect_equal_donor),
                "apa_n_pas": int(apa.n_features),
                "ase_sig": bool(ase.sig),
                "ase_qval": float(ase.qval),
                "ase_effect": float(ase.effect_equal_donor),
                "portal_url": "https://scthread.ai4sc.ac.cn/browse?query=PTPRC",
                "api_path": "/api/gene/PTPRC/overview",
            }
        ]
    ).to_csv(PTPRC, sep="\t", index=False)

    inventory = pd.read_csv(INVENTORY, sep="\t", dtype=str, keep_default_na=False)
    inventory = inventory.set_index("claim_id")
    for axis in ("ase", "diu", "apa"):
        frame = frames[axis]
        upper = axis.upper()
        source = str(P0 / f"{axis}_observed_9999.tsv")
        inventory.loc[f"C1_{axis}_sig", "number"] = str(
            int(frame.sig.astype(bool).sum())
        )
        inventory.loc[f"C1_{axis}_sig", "unit"] = f"genes (of {len(frame)})"
        inventory.loc[f"C1_{axis}_sig", "source_path"] = source
        inventory.loc[f"C1_{axis}_p05", "number"] = (
            f"{float(frame.pval.lt(0.05).mean()):.5f}"
        )
        inventory.loc[f"C1_{axis}_p05", "source_path"] = source
        target = diu if axis == "diu" else apa if axis == "apa" else ase
        detail = (
            f"q={float(target.qval):.6g}; "
            f"equal-donor effect={float(target.effect_equal_donor):.3f}; "
        )
        if axis == "diu":
            detail += f"n_iso={int(target.n_features)}; "
        elif axis == "apa":
            detail += f"n_pas={int(target.n_features)}; "
        detail += f"sig={bool(target.sig)}"
        inventory.loc[f"C3_ptprc_{axis}", "claim_text"] = f"PTPRC {upper}"
        inventory.loc[f"C3_ptprc_{axis}", "number"] = detail
        inventory.loc[f"C3_ptprc_{axis}", "source_path"] = (
            f"{source} gene={GID if axis != 'ase' else 'PTPRC'}"
        )
    inventory.reset_index().to_csv(INVENTORY, sep="\t", index=False)

    print(pd.read_csv(SUMMARY, sep="\t").to_string(index=False))
    print(pd.read_csv(PTPRC, sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
