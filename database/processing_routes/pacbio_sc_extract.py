#!/usr/bin/env python3
# PacBio 10x-3prime single-cell barcode extraction (bypasses lima; works on CCS FASTQ).
# For each read: find 10x R1 adapter (fwd first, then rev-comp), take next 16bp=BC, next 12bp=UMI.
# Keep reads whose BC is on the 10x whitelist. Write FASTQ with CB/UB embedded in read name.
#
# v2 (2026-09-12): HIT-scISOseq / CellRanger-style 1-mismatch cell barcode correction, on by default.
#   Method (docs/processing_logic_20260912/CB_CORRECTION_METHOD.md; Shi et al. 2023 Nat Commun 14:2385,
#   scISA-Tools cellBC_UMI_corrector -> cellranger-cs stats.py correct_bc_error):
#   - prior(w)  = exact-hit count of whitelist barcode w in this run, normalized (pass 1);
#   - candidates = whitelist barcodes at Hamming distance 1 (single substitution) from the
#     observed CB that were observed exactly at least once;
#   - likelihood(w) = prior(w) * 10**(-min(Q_diff,33)/10), Q_diff = base Q at the differing position
#     (Phred+33, capped at 33 as in Cell Ranger);
#   - posterior normalized over candidates; correct to the argmax if posterior > threshold (0.975),
#     else leave unassigned (not written, not a new cell). Unique candidate => posterior 1 => corrected
#     regardless of Q (Cell Ranger behavior; Q only matters when >=2 candidates compete).
#   UMI is NOT corrected (deferred; same as published round 1 scope).
#   --exact-only reproduces v1 (exact-match) behavior byte for byte.
#
# usage: pacbio_sc_extract.py <in.fastq.gz> <whitelist.txt.gz> <out.fastq> <stats.txt>
#                              [--exact-only] [--cb-posterior-threshold 0.975]
import gzip, sys

BC, UMI = 16, 12
adapter = "CTACACGACGCTCTTCCGATCT"; alen = len(adapter)
NUCS = "ACGT"
comp = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(comp)[::-1]

def parse_args(argv):
    pos, opts = [], {"exact_only": False, "threshold": 0.975}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--exact-only": opts["exact_only"] = True
        elif a == "--correct-cb": pass  # default on; explicit no-op
        elif a == "--cb-posterior-threshold":
            i += 1; opts["threshold"] = float(argv[i])
        elif a.startswith("--cb-posterior-threshold="):
            opts["threshold"] = float(a.split("=", 1)[1])
        else: pos.append(a)
        i += 1
    return pos, opts

inf, wlf, outf, statf = parse_args(sys.argv[1:])[0][:4]
opts = parse_args(sys.argv[1:])[1]
THRESH = opts["threshold"]
CORRECT = not opts["exact_only"]

op = gzip.open if wlf.endswith(".gz") else open
wl = set()
with op(wlf, "rt") as w:
    for l in w:
        wl.add(l.strip().split("-")[0])
sys.stderr.write(f"whitelist {len(wl)} barcodes\n")

def modules(seq, qual):
    """Valid adapter-anchored modules on both strands, F first then R (v1 order).
    v1 semantics (must be preserved): a read is kept when EITHER strand's CB exact-hits
    the whitelist — the first such strand wins; F having an adapter does not preclude
    trying R. Correction (v2) extends this: exact match on any strand first, else the
    first strand (F then R) whose CB passes posterior correction."""
    mods = []
    p = seq.find(adapter)
    if p >= 0 and p + alen + BC + UMI <= len(seq):
        mods.append((seq[p+alen:p+alen+BC], qual[p+alen:p+alen+BC],
                     seq[p+alen+BC:p+alen+BC+UMI], seq, qual, "F"))
    r = rc(seq)
    p = r.find(adapter)
    if p >= 0 and p + alen + BC + UMI <= len(r):
        rq = qual[::-1]
        mods.append((r[p+alen:p+alen+BC], rq[p+alen:p+alen+BC],
                     r[p+alen+BC:p+alen+BC+UMI], r, rq, "R"))
    return mods

def correct_bc(cb, cbq, prior, thresh):
    """Cell Ranger correct_bc_error: enumerate Hamming-1 substitutions, keep those on the
    observed-prior whitelist, weight by prior*error-prob at the differing base, normalize,
    return argmax if posterior > thresh. Ties resolve to first candidate in enumeration order."""
    total = 0.0
    best_l = -1.0; best = None; ncand = 0
    for pos in range(len(cb)):
        e = 10.0 ** (-min(ord(cbq[pos]) - 33, 33.0) / 10.0)
        orig = cb[pos]
        for c in NUCS:
            if c == orig: continue
            test = cb[:pos] + c + cb[pos+1:]
            p_bc = prior.get(test)
            if p_bc is not None:
                l = p_bc * e
                total += l; ncand += 1
                if l > best_l: best_l = l; best = test
    if ncand == 0: return None, "unassigned"
    if best_l / total > thresh: return best, "corrected"
    return None, "ambiguous"

def reads(fh):
    while True:
        h = fh.readline()
        if not h: return
        seq = fh.readline().strip(); fh.readline(); qual = fh.readline().strip()
        yield h[1:].split()[0], seq, qual

n = kept = fwd = rev = 0
exact_hits = corrected = ambiguous = unassigned = 0

if CORRECT:
    # pass 1: prior = exact-hit CB counts in this run (first whitelist-hit strand wins, v1 order)
    from collections import Counter
    wl_counts = Counter()
    with gzip.open(inf, "rt") as fh:
        for name, seq, qual in reads(fh):
            exact = next((m for m in modules(seq, qual) if m[0] in wl), None)
            if exact is not None:
                wl_counts[exact[0]] += 1
    wl_sum = sum(wl_counts.values())
    prior = {b: c / wl_sum for b, c in wl_counts.items()} if wl_sum else {}
    sys.stderr.write(f"prior: {len(prior)} observed whitelist barcodes, {wl_sum} exact hits\n")

with gzip.open(inf, "rt") as fh, open(outf, "w") as out:
    for name, seq, qual in reads(fh):
        n += 1
        mods = modules(seq, qual)
        if not mods: continue
        m = next((mm for mm in mods if mm[0] in wl), None)
        if m is not None:
            status = "exact"
        elif CORRECT:
            m = None; any_cand = False
            for mm in mods:
                bc2, st = correct_bc(mm[0], mm[1], prior, THRESH)
                if st == "corrected":
                    mm = (bc2,) + mm[1:]; m = mm; status = "corrected"; break
                if st == "ambiguous": any_cand = True
            if m is None:
                status = "ambiguous" if any_cand else "unassigned"
        else:
            status = "unassigned"
        if status == "exact": exact_hits += 1
        elif status == "corrected": corrected += 1
        elif status == "ambiguous": ambiguous += 1
        else: unassigned += 1
        if status in ("exact", "corrected"):
            bc, bcq, ub, s, q, tag = m
            kept += 1; fwd += (tag == "F"); rev += (tag == "R")
            out.write(f"@{name}|CB:{bc}|UB:{ub}\n{s}\n+\n{q}\n")

mode = "exact-only" if not CORRECT else f"correct-cb"
with open(statf, "w") as st:
    st.write(f"total={n} kept={kept} ({100*kept/max(n,1):.1f}%) fwd={fwd} rev={rev} "
             f"exact_hits={exact_hits} corrected={corrected} ambiguous={ambiguous} "
             f"unassigned={unassigned} cb_correction={mode} threshold={THRESH}\n")
sys.stderr.write(f"total={n} kept={kept} ({100*kept/max(n,1):.1f}%) exact={exact_hits} "
                 f"corrected={corrected} ambiguous={ambiguous} unassigned={unassigned}\n")
