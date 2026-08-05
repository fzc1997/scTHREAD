#!/usr/bin/env python3
"""Mutation tests for scripts/check_anchor_numbers.py.

The control asserts the current prose passes; each mutation reintroduces a
defect that actually reached the repository and asserts the checker rejects it.

The first two cases are the 2026-08-05 incident: the Supplementary Figure 4
legend claimed the published mouse cohort sat outside the frozen release while
the Results said its 68,417 cells were part of the 923,389. Every callout
resolved and all four QA gates passed, because nothing compared what the prose
asserted against the release manifest.

Mutations are guarded against being no-ops, so a reworded manuscript can never
turn a test into a silent pass.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import shutil
import tempfile
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
CHECKER = PROJECT / "scripts/check_anchor_numbers.py"
LEGENDS = PROJECT / "docs/NAR_SF_LEGENDS.md"


def load_checker():
    spec = importlib.util.spec_from_file_location("check_anchor_numbers", CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_case(module, target: Path, label: str, text: str, baseline: str,
             expect_fail: bool) -> None:
    if expect_fail and text == baseline:
        raise AssertionError(
            f"{label}: mutation was a no-op, so this case proves nothing. "
            "The legends have probably been reworded; update the search string.")
    target.write_text(text, encoding="utf-8")
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        code = module.main()
    if bool(code) != expect_fail:
        raise AssertionError(
            f"{label}: expected {'rejection' if expect_fail else 'acceptance'}, "
            f"got exit {code}\n{buffer.getvalue()}")
    print(f"  ok  {label}")


def main() -> None:
    module = load_checker()
    baseline = LEGENDS.read_text(encoding="utf-8")

    sandbox = Path(tempfile.mkdtemp(prefix="anchor_mutation_"))
    try:
        target = sandbox / "legends.md"
        module.PROSE = [module.PROSE[0], target]

        run_case(module, target, "current legends are accepted",
                 baseline, baseline, expect_fail=False)

        run_case(module, target,
                 "claiming the mouse cohort is outside the release is rejected",
                 baseline.replace("The cohort is part of the frozen release:",
                                  "This dataset is outside the frozen release:"),
                 baseline, expect_fail=True)

        run_case(module, target,
                 "claiming its cells are additional to the main text is rejected",
                 baseline.replace(
                     "its 68,417 registered cells are among the 923,389.",
                     "its cells are additional to the counts in the main text."),
                 baseline, expect_fail=True)

        run_case(module, target,
                 "a denominator from a superseded release is rejected",
                 baseline.replace("923,389** cells", "850,938** cells", 1),
                 baseline, expect_fail=True)
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)

    print("ANCHOR-NUMBER MUTATION TESTS PASS")


if __name__ == "__main__":
    main()
