#!/usr/bin/env python3
"""Check the manuscript's anchor numbers against values recomputed from source.

The cross-reference checks prove that every callout resolves. They cannot prove
that what a legend *says* is true. On 2026-08-05 the Supplementary Figure 4
legend and figure footer asserted that the published mouse cohort was outside
the frozen release, while the Results said its 68,417 cells were part of the
923,389. Both files were internally consistent, every callout resolved, and all
four QA gates passed. The registry settles it: those six runs are among the 453
and their cells are among the 923,389, so the legend was wrong.

This script closes that gap for the small set of numbers the whole manuscript
hangs on. Each anchor is recomputed from the release manifest or the served
portal tables and compared against every occurrence in the prose. It also
rejects denominators from superseded releases, which is how the wrong legend
survived a scope change in the first place.

Run standalone, or via manuscript/build_manuscript_docx.slurm.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "figures/scripts"))

PROSE = [
    ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md",
    ROOT / "docs/NAR_SF_LEGENDS.md",
]
SF10_MANIFEST = ROOT / "tables/sf10_v3/sf10_v3_manifest.json"
MOUSE_DATASET = "OWN_ASE_scONT"

# Denominators from superseded releases. Any of these appearing as a live claim
# means a file was written against an older scope and never revisited.
DEAD_VALUES = {
    "434": "run count of the pre-20260803 release",
    "469": "run count of the sample_registry.tsv.bak_2338 snapshot",
    "850,938": "isoquant_cells of the 434-run release",
    "845,781": "isoquant_cells of the bak_2338 snapshot",
    "845,746": "an intermediate cell total",
    "941,733": "a superseded cell total",
    "1,045,231": "run-level sum of isoquant_cells, which double counts studies",
    "452,213": "the labelled analysis subset, which is not a release denominator",
}

# Claims that were true before the scope change and are now false. Targeted
# regression guards, not a general semantic check.
FORBIDDEN_CLAIMS = [
    (r"outside the frozen (?:manuscript )?(?:release|snapshot)",
     f"{MOUSE_DATASET} / CRA044500 is inside the frozen release"),
    (r"(?:its )?cells are additional to the counts in the main text",
     f"{MOUSE_DATASET} / CRA044500 contributes to the main-text counts"),
    (r"is \*?not\*? in the\s+\d+-record manifest",
     f"{MOUSE_DATASET} / CRA044500 is in the release manifest"),
]


def truth() -> dict[str, int]:
    """Recompute every anchor from the data the release is built from."""
    import render_nar_bio as R

    head = R.catalog_headline()
    done = R.registry_done()
    mouse = done[done["gse"].astype(str).eq(MOUSE_DATASET)]
    manifest = json.loads(SF10_MANIFEST.read_text())

    values = {
        "n_runs": int(head["n_runs"]),
        "n_studies": int(head["n_studies"]),
        "n_cells": int(head["n_cells"]),
        "cells_human": int(head["cells_human"]),
        "cells_mouse": int(head["cells_mouse"]),
        "n_ont": int(head["n_ont"]),
        "n_pacbio": int(head["n_pacbio"]),
        "mouse_cohort_runs": int(len(mouse)),
        "mouse_cohort_cells": int(mouse["isoquant_cells"].astype(float).sum()),
        "portal_embedding_cells": int(manifest["portal_embedding"]["cells"]),
        "portal_cell_types": int(manifest["portal_embedding"]["cell_types"]),
    }
    if len(done) != values["n_runs"]:
        raise SystemExit(f"manifest rows {len(done)} != headline runs {values['n_runs']}")
    return values


def main() -> int:
    values = truth()
    problems: list[str] = []

    # The mouse cohort is inside the release, so its cells must be a strict
    # subset of the total and its embedding a subset of its own cells.
    if not values["mouse_cohort_cells"] < values["n_cells"]:
        problems.append("mouse cohort cells are not a subset of the release total")
    if not values["portal_embedding_cells"] <= values["mouse_cohort_cells"]:
        problems.append(
            f"portal embedding ({values['portal_embedding_cells']:,}) exceeds the "
            f"cohort's registered cells ({values['mouse_cohort_cells']:,}); these are "
            "different views of one cohort and must not be swapped")
    if values["cells_human"] + values["cells_mouse"] != values["n_cells"]:
        problems.append("human + mouse cells do not sum to the release total")

    for path in PROSE:
        if not path.is_file():
            problems.append(f"missing prose file: {path}")
            continue
        text = path.read_text(encoding="utf-8")
        body = text.split("## References", 1)[0]
        name = path.name

        for dead, why in DEAD_VALUES.items():
            for m in re.finditer(rf"(?<![\d,]){re.escape(dead)}(?![\d,])", body):
                line = body[:m.start()].count("\n") + 1
                context = body.splitlines()[line - 1].strip()
                # A file may name a dead value in order to forbid it.
                if re.search(r"dead|superseded|stale|forbidden|do not|must not|no longer",
                             context, re.I):
                    continue
                problems.append(f"{name}:{line} uses the dead denominator {dead} "
                                f"({why}): {context[:90]}")

        for pattern, why in FORBIDDEN_CLAIMS:
            for m in re.finditer(pattern, body, re.I):
                line = body[:m.start()].count("\n") + 1
                problems.append(f"{name}:{line} asserts '{m.group(0)}' but {why}")

    print("anchors recomputed from source:")
    for key in sorted(values):
        print(f"  {key:24s} {values[key]:>9,}")
    print(f"\n  mouse cohort {values['mouse_cohort_cells']:,} cells "
          f"({values['mouse_cohort_runs']} runs) is inside the release of "
          f"{values['n_cells']:,} cells ({values['n_runs']} runs)")
    print(f"  its portal embedding holds {values['portal_embedding_cells']:,} "
          f"of those cells")

    if problems:
        print("\nPROBLEMS")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("\nANCHOR NUMBERS AGREE WITH THE RECOMPUTED RELEASE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
