#!/usr/bin/env python3
"""Mutation tests for scripts/check_crossrefs.py.

A cross-reference checker that never fails is indistinguishable from one that
does not run. Each case below breaks the draft in one specific way and asserts
the checker rejects it; the control asserts the unmodified draft still passes.

Every mutation is guarded against being a no-op. That guard exists because the
first version of these tests silently passed: the manuscript had been reworded
underneath a hard-coded search string, so the "mutation" changed nothing and the
checker was credited with catching a defect that was never introduced.

The three failure modes covered are the ones that have actually occurred in this
manuscript:

  * a callout to a supplementary figure that no longer exists (a paragraph
    rewrite on 2026-08-04 reinstated "Supplementary Figures S6 and S7" after the
    pack had been compacted to four figures);
  * supplementary figures cited out of first-mention order, which NAR requires
    and which no earlier check enforced;
  * a panel letter that the figure's legend does not declare.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import shutil
import tempfile
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
CHECKER = PROJECT / "scripts/check_crossrefs.py"


def load_checker():
    spec = importlib.util.spec_from_file_location("check_crossrefs", CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_case(module, sandbox: Path, label: str, text: str, baseline: str,
             expect_fail: bool) -> None:
    if expect_fail and text == baseline:
        raise AssertionError(
            f"{label}: mutation was a no-op, so this case proves nothing. "
            "The draft has probably been reworded; update the search string.")
    draft = sandbox / "draft.md"
    draft.write_text(text, encoding="utf-8")
    module.DRAFT = draft
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        code = module.main()
    failed = bool(code)
    if failed != expect_fail:
        raise AssertionError(
            f"{label}: expected {'rejection' if expect_fail else 'acceptance'}, "
            f"got exit {code}\n{buffer.getvalue()}")
    print(f"  ok  {label}")


def main() -> None:
    module = load_checker()
    baseline = Path(module.DRAFT).read_text(encoding="utf-8")

    sandbox = Path(tempfile.mkdtemp(prefix="crossref_mutation_"))
    try:
        figdir = sandbox / "figures"
        figdir.mkdir()
        for n in (1, 2, 3, 4):
            (figdir / f"NAR_SF{n}.pdf").write_bytes(b"")
        module.FIGDIR = figdir

        run_case(module, sandbox, "unmodified draft is accepted",
                 baseline, baseline, expect_fail=False)

        # Move the first supplementary-figure callout to S3 so that S3 is
        # first-mentioned ahead of S1.
        head = baseline.index("Supplementary Figure S1")
        span = len("Supplementary Figure S1")
        run_case(module, sandbox, "out-of-order first mention is rejected",
                 baseline[:head] + "Supplementary Figure S3" + baseline[head + span:],
                 baseline, expect_fail=True)

        run_case(module, sandbox, "callout to a deleted figure is rejected",
                 baseline.replace("Supplementary Figure S2e",
                                  "Supplementary Figures S6 and S7"),
                 baseline, expect_fail=True)

        # Supplementary Figure 2 declares panels a-e in its legend.
        panel = baseline.index("Supplementary Figure S2c")
        span = len("Supplementary Figure S2c")
        run_case(module, sandbox, "undeclared panel letter is rejected",
                 baseline[:panel] + "Supplementary Figure S2f" + baseline[panel + span:],
                 baseline, expect_fail=True)

        # An artefact nothing points at is a leftover from a renumbering.
        (figdir / "NAR_SF7.pdf").write_bytes(b"")
        run_case(module, sandbox, "uncited figure left in figures/ is rejected",
                 baseline + "\n", baseline, expect_fail=True)
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)

    print("CROSS-REFERENCE MUTATION TESTS PASS")


if __name__ == "__main__":
    main()
