#!/usr/bin/env python3
"""
F2.1 feature component — splice-site strength via position-weight matrices (self-contained).

Builds 5'SS (donor) and 3'SS (acceptor) PWMs from the reference junction set (1.36M GENCODE
junctions) + genome FASTA, then scores any junction's donor/acceptor log-odds strength.
First-principles (no MaxEntScan / external tool). Windows follow the standard MaxEnt convention:
  5'SS donor:    exon[-3..-1] + intron[+1..+6]   (9 nt around the exon|intron boundary)
  3'SS acceptor: intron[-20..-1] + exon[+1..+3]   (23 nt around the intron|exon boundary)

Coordinate convention (matches F2 chain elems "A-B", A=exon_i_end, B=exon_{i+1}_start, 1-based):
  intron 1-based = [A+1, B-1]; genome 0-based slice.
  donor (+ strand)   = 9nt at exon|intron boundary around A
  acceptor (+ strand)= 23nt at intron|exon boundary around B
  - strand: reverse-complement, donor/acceptor swap.

Usage: build PWM once -> save; import score_donor(chrom,A,strand)/score_acceptor(chrom,B,strand).
Run standalone to BUILD + UNIT-TEST: f2_splice_pwm.py build
"""
import sys, os, json
import numpy as np

FASTA = os.environ.get("SCTHREAD_GENOME_FASTA", "reference/genome.fa")
REFJCT = os.environ.get("SCTHREAD_REFERENCE_JUNCTIONS", "reference/ref_genes.gencode_v44_filtered.junctions.bed")
PWM_OUT = os.environ.get("SCTHREAD_SPLICE_PWM", "outputs/f2_grammar/splice_pwm.json")

BASES = "ACGT"; BIDX = {b: i for i, b in enumerate(BASES)}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
def rc(s): return "".join(COMP.get(c, "N") for c in reversed(s))

# donor window: 3 exonic + 6 intronic = 9 nt; acceptor: 20 intronic + 3 exonic = 23 nt
DON_EX, DON_IN = 3, 6
ACC_IN, ACC_EX = 20, 3

def _fa():
    from pyfaidx import Fasta
    return Fasta(FASTA, rebuild=False)

def donor_seq(fa, chrom, A, strand):
    # BED/ref intron: start(0-based)=A, so donor at exon|intron boundary = genome[A-3 : A+6] (+ strand)
    if strand == "+":
        s = str(fa[chrom][A - DON_EX: A + DON_IN]).upper()
    else:
        s = rc(str(fa[chrom][A - DON_IN: A + DON_EX]).upper())
    return s

def acceptor_seq(fa, chrom, B, strand):
    # BED intron end(exclusive)=B; acceptor at intron|exon boundary = genome[B-20 : B+3] (+ strand)
    if strand == "+":
        s = str(fa[chrom][B - ACC_IN: B + ACC_EX]).upper()
    else:
        s = rc(str(fa[chrom][B - ACC_EX: B + ACC_IN]).upper())
    return s

def build_pwm(seqs, L):
    cnt = np.ones((L, 4))  # pseudocount 1
    n = 0
    for s in seqs:
        if len(s) != L or "N" in s:
            continue
        for i, b in enumerate(s):
            cnt[i, BIDX[b]] += 1
        n += 1
    freq = cnt / cnt.sum(axis=1, keepdims=True)
    return freq, n

def score(seq, freq, bg=0.25):
    if len(seq) != freq.shape[0] or "N" in seq:
        return float("nan")
    return float(sum(np.log2(freq[i, BIDX[b]] / bg) for i, b in enumerate(seq)))

def build():
    import duckdb
    fa = _fa()
    # sample reference junctions to build PWM (use up to 200k for speed; BED: chrom,start,end,...,strand col6)
    don_seqs, acc_seqs = [], []
    with open(REFJCT) as fh:
        for k, line in enumerate(fh):
            if k >= 200000:
                break
            f = line.split("\t")
            if len(f) < 6:
                continue
            chrom, start, end, strand = f[0], int(f[1]), int(f[2]), f[5].strip()
            if chrom not in fa:
                continue
            d = donor_seq(fa, chrom, start, strand)      # BED start = intron 0-based start = donor anchor A
            a = acceptor_seq(fa, chrom, end, strand)      # BED end = intron exclusive end = acceptor anchor B
            if d: don_seqs.append(d)
            if a: acc_seqs.append(a)
    don_freq, nd = build_pwm(don_seqs, DON_EX + DON_IN)
    acc_freq, na = build_pwm(acc_seqs, ACC_IN + ACC_EX)
    json.dump({"donor": don_freq.tolist(), "acceptor": acc_freq.tolist(),
               "don_L": DON_EX + DON_IN, "acc_L": ACC_IN + ACC_EX, "n_donor": nd, "n_acceptor": na,
               "windows": {"DON_EX": DON_EX, "DON_IN": DON_IN, "ACC_IN": ACC_IN, "ACC_EX": ACC_EX}},
              open(PWM_OUT, "w"))
    print(f"built PWM from {nd:,} donor / {na:,} acceptor reference splice sites -> {PWM_OUT}")
    return don_freq, acc_freq

def unit_test():
    d = json.load(open(PWM_OUT))
    don_freq = np.array(d["donor"]); acc_freq = np.array(d["acceptor"])
    fa = _fa()
    # canonical GT-AG reference junctions should score HIGH; shuffled/random positions LOW
    import random
    real_d, real_a, rand_d = [], [], []
    with open(REFJCT) as fh:
        lines = [next(fh) for _ in range(5000)]
    for line in lines[:2000]:
        f = line.split("\t")
        if len(f) < 6: continue
        chrom, start, end, strand = f[0], int(f[1]), int(f[2]), f[5].strip()
        if chrom not in fa: continue
        real_d.append(score(donor_seq(fa, chrom, start, strand), don_freq))
        real_a.append(score(acceptor_seq(fa, chrom, end, strand), acc_freq))
        # random position on same chrom
        rp = start + random.randint(5000, 50000)
        rand_d.append(score(donor_seq(fa, chrom, rp, strand), don_freq))
    real_d = [x for x in real_d if x == x]; rand_d = [x for x in rand_d if x == x]
    real_a = [x for x in real_a if x == x]
    print(f"UNIT TEST: donor score real={np.mean(real_d):.2f} vs random={np.mean(rand_d):.2f} "
          f"(real should be >> random)")
    print(f"           acceptor score real={np.mean(real_a):.2f}")
    ok = np.mean(real_d) > np.mean(rand_d) + 3
    print("  PASS" if ok else "  FAIL — donor PWM not discriminating")
    return ok

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "build":
        os.makedirs(os.path.dirname(PWM_OUT), exist_ok=True)
        build()
        unit_test()
