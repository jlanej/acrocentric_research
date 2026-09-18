"""
Transmission of acrocentric short-arm block structure through the CEPH1463 pedigree.

Asks how good a reference a parent's assembled short arm is for their child's, using
only the satellite block annotation. Because the arms do not cross over, a child's arm
should be a copy of one parental arm; the measurement is whether that shows up as
block-level identity, and how it compares with the non-transmitting parent's arms.

Design. All arms are anchored at the proximal (centromeric) end, which is the end that
assembles; the distal end is truncated by variable amounts. Comparisons therefore use a
semi-global alignment with free terminal gaps at the distal end only. Two metrics are
reported per pair: order identity (matched columns / aligned columns) and size
divergence (median |log2 length ratio| over matched columns).

Matched-candidate contrast. For each child arm the transmitting parent is compared with
the non-transmitting parent, each contributing at most two candidate arms, so the
minimum is taken over the same number of candidates on both sides. NA12877 and NA12878
are unrelated, which makes the non-transmitting parent an internal population control.

Usage:  python pedigree_transmission.py [datadir] [outdir]
"""
import csv
import json
import os
import statistics
import sys
from collections import defaultdict

MISMATCH = 1.0
GAP = 1.0
MIN_BLOCKS = 5
MIN_MATCH_FOR_SIZE = 3

FATHER, MOTHER = "NA12877", "NA12878"
G3 = ["NA12879", "NA12881", "NA12882", "NA12883",
      "NA12884", "NA12885", "NA12886", "NA12887"]
CHROMS = ["chr13", "chr14", "chr15", "chr21", "chr22"]


# ---------------------------------------------------------------- loading

def load_arms(datadir):
    """Return {(sample, hap, chrom): [(family, length_bp), ...]} proximal-first."""
    meta = {}
    for r in csv.DictReader(open(os.path.join(datadir, "arm_syntax_pedigree.csv"))):
        meta[r["tig"]] = (r["sample"], int(r["hap"]), r["chrom"], r["syntax"])
    raw = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(datadir, "arm_blocks_pedigree.csv"))):
        raw[r["tig"]].append((int(r["start"]), int(r["end"]), r["fam"]))
    arms, flipped = {}, 0
    for tig, blocks in raw.items():
        if tig not in meta:
            continue
        sample, hap, chrom, syntax = meta[tig]
        blocks.sort()
        fams = [b[2] for b in blocks]
        # the syntax field is written telomere -> centromere; detect tig orientation
        if "-".join(fams) != syntax:
            if "-".join(reversed(fams)) == syntax:
                blocks = blocks[::-1]
                flipped += 1
            else:
                continue  # block set and syntax field disagree; skip
        # store proximal-first: reverse the telomere -> centromere order
        arms[(sample, hap, chrom)] = [(b[2], b[1] - b[0]) for b in blocks][::-1]
    return arms, flipped


# ---------------------------------------------------------------- alignment

def align_proximal(a, b):
    """Semi-global alignment, gaps free at the distal (tail) end only.

    a, b are proximal-first [(family, length)] lists. Returns
    (order_identity, size_divergence, n_cols, n_match) or None if uninformative.
    """
    n, m = len(a), len(b)
    if n < MIN_BLOCKS or m < MIN_BLOCKS:
        return None
    INF = float("inf")
    s = [[INF] * (m + 1) for _ in range(n + 1)]
    bt = [[0] * (m + 1) for _ in range(n + 1)]
    s[0][0] = 0.0
    for i in range(1, n + 1):
        s[i][0] = i * GAP
        bt[i][0] = 1
    for j in range(1, m + 1):
        s[0][j] = j * GAP
        bt[0][j] = 2
    for i in range(1, n + 1):
        ai = a[i - 1][0]
        si, sp = s[i], s[i - 1]
        bi = bt[i]
        for j in range(1, m + 1):
            d = sp[j - 1] + (0.0 if ai == b[j - 1][0] else MISMATCH)
            u = sp[j] + GAP
            l = si[j - 1] + GAP
            if d <= u and d <= l:
                si[j], bi[j] = d, 0
            elif u <= l:
                si[j], bi[j] = u, 1
            else:
                si[j], bi[j] = l, 2
    # free terminal gaps: best cell on the last row or last column
    best, bi_, bj_ = INF, n, m
    for j in range(MIN_BLOCKS, m + 1):
        if s[n][j] < best:
            best, bi_, bj_ = s[n][j], n, j
    for i in range(MIN_BLOCKS, n + 1):
        if s[i][m] < best:
            best, bi_, bj_ = s[i][m], i, m
    i, j = bi_, bj_
    cols = matches = 0
    ratios = []
    while i > 0 and j > 0:
        t = bt[i][j]
        if t == 0:
            cols += 1
            if a[i - 1][0] == b[j - 1][0]:
                matches += 1
                la, lb = a[i - 1][1], b[j - 1][1]
                if la > 0 and lb > 0:
                    ratios.append(abs(_log2(la / lb)))
            i, j = i - 1, j - 1
        elif t == 1:
            cols += 1
            i -= 1
        else:
            cols += 1
            j -= 1
    if cols < MIN_BLOCKS:
        return None
    ident = matches / cols
    size_div = statistics.median(ratios) if len(ratios) >= MIN_MATCH_FOR_SIZE else None
    return ident, size_div, cols, matches


def _log2(x):
    import math
    return math.log2(x)


def best_against(arm, cands):
    """Best (highest order identity, then lowest size divergence) over candidates."""
    best = None
    for c in cands:
        r = align_proximal(arm, c)
        if r is None:
            continue
        key = (r[0], -(r[1] if r[1] is not None else 9.9))
        if best is None or key > best[0]:
            best = (key, r)
    return None if best is None else best[1]


# ---------------------------------------------------------------- analyses

def transmission_contrast(arms):
    """Child arms vs transmitting and non-transmitting parent, matched candidates."""
    rows = []
    for child in G3:
        for chrom in CHROMS:
            fa = [arms[(FATHER, h, chrom)] for h in (1, 2) if (FATHER, h, chrom) in arms]
            mo = [arms[(MOTHER, h, chrom)] for h in (1, 2) if (MOTHER, h, chrom) in arms]
            kids = [(h, arms[(child, h, chrom)]) for h in (1, 2)
                    if (child, h, chrom) in arms]
            if not fa or not mo or not kids:
                continue
            for hap, arm in kids:
                rf, rm = best_against(arm, fa), best_against(arm, mo)
                if rf is None or rm is None:
                    continue
                pat = rf[0] > rm[0] or (rf[0] == rm[0] and (rf[1] or 9) < (rm[1] or 9))
                rows.append(dict(
                    child=child, chrom=chrom, hap=hap, n_blocks=len(arm),
                    assigned="paternal" if pat else "maternal",
                    ident_transmitting=round(rf[0] if pat else rm[0], 4),
                    ident_other=round(rm[0] if pat else rf[0], 4),
                    size_div_transmitting=_r(rf[1] if pat else rm[1]),
                    size_div_other=_r(rm[1] if pat else rf[1]),
                    cols=rf[2] if pat else rm[2]))
    return rows


def parent_of_origin_consistency(rows):
    """Do a child's two arms of a chromosome resolve to one paternal + one maternal?"""
    by = defaultdict(list)
    for r in rows:
        by[(r["child"], r["chrom"])].append(r["assigned"])
    pairs = {k: v for k, v in by.items() if len(v) == 2}
    ok = sum(1 for v in pairs.values() if set(v) == {"paternal", "maternal"})
    return dict(chromosomes_with_both_haplotypes=len(pairs), one_of_each=ok,
                fraction=round(ok / len(pairs), 4) if pairs else None)


def population_order_baseline(datadir):
    """Order identity between unrelated HPRC arms of the same chromosome (context)."""
    arms = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(datadir, "arm_inventory.csv"))):
        if r["source"] != "censat":
            continue
        fams = r["syntax"].split("-")
        if len(fams) >= MIN_BLOCKS:
            arms[r["chrom"]].append((r["sample"], [(f, 1) for f in fams][::-1]))
    import random
    rng = random.Random(7)
    out = []
    for chrom, lst in arms.items():
        rng.shuffle(lst)
        lst = lst[:60]
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                if lst[i][0] == lst[j][0]:
                    continue
                r = align_proximal(lst[i][1], lst[j][1])
                if r:
                    out.append(dict(chrom=chrom, ident=round(r[0], 4)))
    return out


def g4_grandparental(arms):
    """Fourth-generation arms vs the G2 couple, and vs each G3 individual.

    Each G4 individual has one married-in parent absent from the data, so about half of
    their arms should have no match in the pedigree: an internal negative control.
    """
    g4 = sorted({s for (s, _, _) in arms if s.startswith("K2000")})
    rows, parent_votes = [], defaultdict(lambda: defaultdict(int))
    for ind in g4:
        for chrom in CHROMS:
            g2 = [arms[(p, h, chrom)] for p in (FATHER, MOTHER) for h in (1, 2)
                  if (p, h, chrom) in arms]
            for hap in (1, 2):
                arm = arms.get((ind, hap, chrom))
                if arm is None or not g2:
                    continue
                r2 = best_against(arm, g2)
                if r2 is None:
                    continue
                scored = {}
                for g3 in G3:
                    c = [arms[(g3, h, chrom)] for h in (1, 2) if (g3, h, chrom) in arms]
                    if not c:
                        continue
                    r3 = best_against(arm, c)
                    if r3:
                        scored[g3] = (r3[0], r3[1] if r3[1] is not None else 9.9)
                tied, best_id = [], -1.0
                if scored:
                    best_id = max(v[0] for v in scored.values())
                    at_top = {g: v for g, v in scored.items() if v[0] == best_id}
                    best_sd = min(v[1] for v in at_top.values())
                    tied = sorted(g for g, v in at_top.items() if v[1] <= best_sd + 0.05)
                rows.append(dict(individual=ind, chrom=chrom, hap=hap,
                                 ident_vs_G2=round(r2[0], 4),
                                 size_div_vs_G2=_r(r2[1]),
                                 ident_vs_best_G3=round(best_id, 4),
                                 n_G3_tied=len(tied),
                                 G3_tied=";".join(tied)))
                if tied and best_id >= 0.9:
                    for g in tied:
                        parent_votes[ind][g] += 1
    votes = []
    n_arms_with_tie = defaultdict(int)
    for r in rows:
        if r["ident_vs_best_G3"] >= 0.9:
            n_arms_with_tie[r["individual"]] += 1
    for ind, v in parent_votes.items():
        nmatched = n_arms_with_tie[ind]
        top = max(v.values())
        winners = sorted(g for g, c in v.items() if c == top)
        votes.append(dict(individual=ind, matched_arms=nmatched,
                          top_vote_count=top, n_top_tied=len(winners),
                          top_candidates=";".join(winners),
                          n_G3_ever_matched=len(v)))
    return rows, sorted(votes, key=lambda d: d["individual"])


def _r(x):
    return None if x is None else round(x, 4)


# ---------------------------------------------------------------- main

def main(datadir, outdir):
    arms, flipped = load_arms(datadir)
    print("[1] pedigree arms loaded: %d (%d tigs reversed to telomere-first)"
          % (len(arms), flipped))

    rows = transmission_contrast(arms)
    _write(os.path.join(outdir, "transmission_contrast.csv"), rows)
    it = [r["ident_transmitting"] for r in rows]
    io = [r["ident_other"] for r in rows]
    st = [r["size_div_transmitting"] for r in rows if r["size_div_transmitting"] is not None]
    so = [r["size_div_other"] for r in rows if r["size_div_other"] is not None]
    print("[2] child arms compared: %d" % len(rows))
    print("    order identity   transmitting %.3f   non-transmitting %.3f"
          % (statistics.median(it), statistics.median(io)))
    print("    exact order      transmitting %.1f%%  non-transmitting %.1f%%"
          % (100 * sum(x == 1.0 for x in it) / len(it),
             100 * sum(x == 1.0 for x in io) / len(io)))
    print("    size divergence  transmitting %.3f   non-transmitting %.3f  (median |log2|)"
          % (statistics.median(st), statistics.median(so)))

    poo = parent_of_origin_consistency(rows)
    print("[3] parent of origin: %d/%d chromosomes give one paternal + one maternal (%s)"
          % (poo["one_of_each"], poo["chromosomes_with_both_haplotypes"], poo["fraction"]))

    g4rows, votes = g4_grandparental(arms)
    _write(os.path.join(outdir, "g4_grandparental.csv"), g4rows)
    _write(os.path.join(outdir, "g4_parent_votes.csv"), votes)
    hi = [r for r in g4rows if r["ident_vs_G2"] >= 0.9]
    print("[4] G4 arms: %d; %d (%.0f%%) match a G2 arm at >=0.9 order identity"
          % (len(g4rows), len(hi), 100 * len(hi) / max(1, len(g4rows))))
    tie = [r["n_G3_tied"] for r in g4rows if r["ident_vs_best_G3"] >= 0.9]
    if tie:
        print("    G3 individuals indistinguishable per matched arm: median %d of 8"
              % statistics.median(tie))
    for v in votes:
        print("    %s: %d matched arms, best G3 vote %d, tied candidates %d (%s)"
              % (v["individual"], v["matched_arms"], v["top_vote_count"],
                 v["n_top_tied"], v["top_candidates"]))

    base = population_order_baseline(datadir)
    _write(os.path.join(outdir, "population_order_baseline.csv"), base)
    print("[5] unrelated HPRC same-chromosome order identity: median %.3f (n=%d pairs)"
          % (statistics.median([r["ident"] for r in base]), len(base)))

    summ = dict(
        n_child_arms=len(rows),
        median_order_identity_transmitting=round(statistics.median(it), 4),
        median_order_identity_non_transmitting=round(statistics.median(io), 4),
        pct_exact_order_transmitting=round(100 * sum(x == 1.0 for x in it) / len(it), 2),
        pct_exact_order_non_transmitting=round(100 * sum(x == 1.0 for x in io) / len(io), 2),
        median_size_divergence_transmitting=round(statistics.median(st), 4),
        median_size_divergence_non_transmitting=round(statistics.median(so), 4),
        parent_of_origin=poo,
        n_g4_arms=len(g4rows),
        pct_g4_matching_G2=round(100 * len(hi) / max(1, len(g4rows)), 2),
        median_unrelated_order_identity=round(
            statistics.median([r["ident"] for r in base]), 4))
    json.dump(summ, open(os.path.join(outdir, "transmission_summary.json"), "w"), indent=1)
    print("\ndone; tables in %s" % outdir)
    return summ


def _write(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "../data"
    o = sys.argv[2] if len(sys.argv) > 2 else d
    main(d, o)
