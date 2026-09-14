#!/usr/bin/env python3
# Kinnex/MAS 10x-3prime single-cell extraction from CCS FASTQ (bypasses skera + lima).
# Each concatenated array read carries N cDNAs (~16 for Kinnex-16), each with its own
# 10x R1 adapter -> 16bp CB -> 12bp UMI -> cDNA. We split the read at EVERY adapter and
# emit one FASTQ record per segment (cDNA), CB/UB in the read name — same format as
# pacbio_sc_extract.py, so the downstream minimap2->tag->IsoQuant flow is unchanged.
#
# v2 (2026-09-12): HIT-scISOseq / CellRanger-style 1-mismatch CB correction per segment,
#   on by default (same algorithm as pacbio_sc_extract.py v2; see
#   docs/processing_logic_20260912/CB_CORRECTION_METHOD.md). Prior = exact-hit CB counts of
#   this run (per segment); candidates = observed whitelist barcodes at Hamming distance 1;
#   likelihood = prior * 10**(-min(Q_diff,33)/10); posterior normalized over candidates;
#   argmax wins if posterior > threshold (0.975), else segment stays unassigned (not written).
#   UMI is NOT corrected (deferred). --exact-only reproduces v1 behavior byte for byte.
#
# usage: pacbio_sc_extract_kinnex.py <in.fastq.gz> <whitelist.txt.gz> <out.fastq> <stats.txt>
#                                    [--exact-only] [--cb-posterior-threshold 0.975]
import gzip, sys

adapter = "CTACACGACGCTCTTCCGATCT"; alen = len(adapter)
BC, UMI = 16, 12
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

pos_args, opts = parse_args(sys.argv[1:])
inf, wlf, outf, statf = pos_args[:4]
THRESH = opts["threshold"]
CORRECT = not opts["exact_only"]

op = gzip.open if wlf.endswith(".gz") else open
wl = set()
with op(wlf, "rt") as w:
    for l in w:
        wl.add(l.strip().split("-")[0])
sys.stderr.write(f"whitelist {len(wl)} barcodes\n")

def all_pos(s):
    out = []; i = s.find(adapter)
    while i >= 0:
        out.append(i); i = s.find(adapter, i + 1)
    return out

def segments(seq, qual):
    """Yield (j, bc, bc_qual, ub, cdna, cdna_qual) per adapter-delimited segment on the
    array strand (more adapter hits wins, ties -> fwd), replicating v1 exactly."""
    fpos = all_pos(seq)
    rseq = rc(seq); rqual = qual[::-1]; rpos = all_pos(rseq)
    if len(rpos) > len(fpos):
        s, q, pos = rseq, rqual, rpos
    else:
        s, q, pos = seq, qual, fpos
    pos.append(len(s))  # sentinel: last segment runs to end
    for j in range(len(pos) - 1):
        p = pos[j]
        if p + alen + BC + UMI > len(s):
            continue
        bc = s[p+alen:p+alen+BC]; ub = s[p+alen+BC:p+alen+BC+UMI]
        cdna_start = p + alen + BC + UMI
        cdna_end = pos[j+1]                 # up to next adapter (or read end)
        cd = s[cdna_start:cdna_end]; cq = q[cdna_start:cdna_end]
        yield j, bc, q[p+alen:p+alen+BC], ub, cd, cq

def correct_bc(cb, cbq, prior, thresh):
    """Cell Ranger correct_bc_error (see pacbio_sc_extract.py v2 docstring)."""
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

n_reads = 0
exact_hits = corrected = ambiguous = unassigned = 0

if CORRECT:
    # pass 1: prior = per-segment exact-hit CB counts in this run
    from collections import Counter
    wl_counts = Counter()
    with gzip.open(inf, "rt") as fh:
        for name, seq, qual in reads(fh):
            for j, bc, bcq, ub, cd, cq in segments(seq, qual):
                if bc in wl:
                    wl_counts[bc] += 1
    wl_sum = sum(wl_counts.values())
    prior = {b: c / wl_sum for b, c in wl_counts.items()} if wl_sum else {}
    sys.stderr.write(f"prior: {len(prior)} observed whitelist barcodes, {wl_sum} exact segment hits\n")

n_seg = kept = 0
seg_per_read = []
with gzip.open(inf, "rt") as fh, open(outf, "w") as out:
    for name, seq, qual in reads(fh):
        n_reads += 1
        k_this = 0
        for j, bc, bcq, ub, cd, cq in segments(seq, qual):
            n_seg += 1
            if len(cd) < 50:            # skip empty/degenerate segments (v1 rule)
                continue
            if bc in wl:
                status = "exact"
            elif CORRECT:
                bc2, status = correct_bc(bc, bcq, prior, THRESH)
                if status == "corrected": bc = bc2
            else:
                status = "unassigned"
            if status == "exact": exact_hits += 1
            elif status == "corrected": corrected += 1
            elif status == "ambiguous": ambiguous += 1
            else: unassigned += 1
            if status not in ("exact", "corrected"):
                continue
            kept += 1; k_this += 1
            out.write(f"@{name}.s{j}|CB:{bc}|UB:{ub}\n{cd}\n+\n{cq}\n")
        seg_per_read.append(k_this)

mean_seg = sum(seg_per_read) / max(len(seg_per_read), 1)
mode = "exact-only" if not CORRECT else "correct-cb"
with open(statf, "w") as st:
    st.write(f"reads={n_reads} segments_kept={kept} mean_seg_per_read={mean_seg:.2f} "
             f"segments_seen={n_seg} exact_hits={exact_hits} corrected={corrected} "
             f"ambiguous={ambiguous} unassigned={unassigned} cb_correction={mode} threshold={THRESH}\n")
sys.stderr.write(f"reads={n_reads} segments_kept={kept} mean_seg_per_read={mean_seg:.2f} "
                 f"exact={exact_hits} corrected={corrected} ambiguous={ambiguous} unassigned={unassigned}\n")
