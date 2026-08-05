#!/usr/bin/env python3
"""Save the current manuscript build under a version name in manuscript/versions/.

The build writes a fixed filename and overwrites it in place, so without this
step the only copies of earlier deliveries live on the iMac. Snapshot first,
then hand the versioned paths straight to the OneDrive sync helper - no staging
directory and no renaming step.

    python3 scripts/snapshot_manuscript_version.py --tag review10_xyz
    python3 .../sync_exact.py --target-dir '<...>' manuscript/versions/<the files>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSIONS = ROOT / "manuscript/versions"
SOURCES = [
    ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.docx",
    ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md",
]
LEDGER = VERSIONS / "VERSIONS.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_version_advances(tag: str) -> None:
    """A Vn tag must be the next number, not one already spent.

    Two numbering schemes ran side by side here - the collaborator's V1/V2/V3
    and the per-edit review tags - and a snapshot once went out as V4 when the
    ledger was already at V8. A reader given a lower number reads it as older
    work. Only Vn tags are policed; review-style tags carry no ordering.
    """
    match = re.match(r"V(\d+)(?![0-9])", tag)
    if not match or not LEDGER.is_file():
        return
    spent = {int(m.group(1))
             for line in LEDGER.read_text().splitlines()[1:]
             if (m := re.match(r"V(\d+)(?![0-9])", line.split("\t", 1)[0]))}
    if not spent:
        return
    wanted, highest = int(match.group(1)), max(spent)
    # re-running the SAME tag adds files to an existing snapshot and is safe:
    # the copy step already refuses to overwrite a file whose content differs.
    # What must be blocked is a version that goes backward, or a new tag that
    # reuses a number already spent under a different name.
    existing = {line.split("\t", 1)[0] for line in LEDGER.read_text().splitlines()[1:]}
    if tag in existing:
        return
    if wanted <= highest:
        raise SystemExit(
            f"tag {tag} reuses or goes back from V{highest}, the highest in "
            f"{LEDGER.name}; the next version is V{highest + 1}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True,
                        help="version tag, e.g. review10_reviewer_response")
    parser.add_argument("--date", default=date.today().isoformat().replace("-", ""))
    parser.add_argument("--note", default="", help="one line recorded in the ledger")
    parser.add_argument("--extra", type=Path, nargs="*", default=[],
                        help="additional files to snapshot under the same tag")
    args = parser.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.tag):
        raise SystemExit("tag must be filename-safe")
    VERSIONS.mkdir(parents=True, exist_ok=True)
    check_version_advances(args.tag)

    written = []
    for src in [*SOURCES, *args.extra]:
        if not src.is_file():
            raise SystemExit(f"missing source: {src}")
        stem = "scTHREAD_NAR_Database_Issue" if "draft" in src.stem else src.stem
        dst = VERSIONS / f"{stem}_{args.tag}_{args.date}{src.suffix}"
        if dst.exists() and sha256(dst) != sha256(src):
            raise SystemExit(f"refusing to overwrite a different file: {dst}")
        shutil.copy2(src, dst)
        written.append(dst)

    rows = []
    if LEDGER.is_file():
        rows = LEDGER.read_text().splitlines()
    if not rows:
        rows = ["tag\tdate\tfile\tbytes\tsha256\tnote"]
    for dst in written:
        rows.append(f"{args.tag}\t{args.date}\t{dst.name}\t{dst.stat().st_size}"
                    f"\t{sha256(dst)}\t{args.note}")
    LEDGER.write_text("\n".join(rows) + "\n")

    print(json.dumps({"tag": args.tag,
                      "files": [str(d.relative_to(ROOT)) for d in written],
                      "ledger": str(LEDGER.relative_to(ROOT))}, indent=2))
    print("\nsync these paths directly - no staging directory, no rename:")
    for d in written:
        print(f"  {d}")


if __name__ == "__main__":
    main()
