#!/usr/bin/env python
"""
Reproduce the satellite-syntax analysis of acrocentric short arms.

Inputs (see docs/DATA.md for how to fetch them):
  <work>/all_samples_censat.bed          Platinum Pedigree per-haplotype satellite track
  <work>/flagger_nucfreq_merged.bed      assembly-error exclusion set for the same arms
  <work>/hprc/*.cenSat.bed               HPRC release 2 per-haplotype CenSat annotations
  <work>/aux/*.fai                       scaffold lengths for those assemblies

Outputs: result tables into <out>/ and figures into <fig>/.

Usage:  python run_analysis.py --work ../work --out ../data --fig ../figures
"""
from __future__ import annotations
import argparse, collections, itertools, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import armsyntax as AS

CHR = ["chr13", "chr14", "chr15", "chr21", "chr22"]
SATS = ["bSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "gSat"]
SATS2 = SATS + ["SST1", "ACRO"]
MASS_MIN = 0.25
SEED = 7


# ---------------------------------------------------------------- utilities
def predict_slot(g, j, alphabet=SATS2, mass_min=MASS_MIN):
    """Family prediction for a block placed in slot j.

    The slot's distribution is renormalised over `alphabet`, but only if the
    slot carries at least `mass_min` total probability on that alphabet.
    Without the gate, a slot dominated by non-satellite sequence would still
    return its highest-probability satellite family from a vanishing
    probability, which is meaningless.
    """
    if j < 0:
        return None, np.nan, 0.0
    pv = np.array([g.E[j, g.si[x]] for x in alphabet])
    mass = float(pv.sum())
    if mass < mass_min:
        return None, np.nan, mass
    pn = pv / mass
    return alphabet[int(np.argmax(pn))], float(pn.max()), mass


def segment_records(segs, ids):
    out = []
    for i in ids:
        a, d, p = segs[i]
        if len(d) >= 4:
            out.append((a, "distal", [b[0] for b in d], d))
        if len(p) >= 4:
            out.append((a, "proximal", [b[0] for b in p], p))
    return out


def fit_models(segs, train_ids, rounds=10):
    tr = segment_records(segs, train_ids)
    M = {}
    dl = [r[2] for r in tr if r[1] == "distal"]
    if len(dl) >= 20:
        M[("distal", "pooled")] = AS.SlotGrammar(seed=AS.seed_from_medoid(dl)).fit(dl, rounds=rounds)
    for c in CHR:
        pl = [r[2] for r in tr if r[1] == "proximal" and r[0].chrom == c]
        if len(pl) >= 20:
            M[("proximal", c)] = AS.SlotGrammar(seed=AS.seed_from_medoid(pl)).fit(pl, rounds=rounds)
    cnt = collections.Counter(); posb = collections.defaultdict(collections.Counter)
    for a, seg, sy, bl in tr:
        n = len(bl)
        for k, (f, s, e) in enumerate(bl):
            if f in SATS:
                cnt[f] += 1
                posb[(seg, min(int(10 * k / max(n, 1)), 9))][f] += 1
    maj = cnt.most_common(1)[0][0]
    posmod = {k: v.most_common(1)[0][0] for k, v in posb.items()}
    return M, maj, posmod


def mask_and_recover(segs, M, maj, posmod, test_ids, rng, alphabet=SATS, max_arms=400):
    """Mask one block at a time and try to recover its family from context."""
    te = segment_records(segs, test_ids); rng.shuffle(te); te = te[:max_arms]
    rows = []
    for a, seg, sy, bl in te:
        key = ("distal", "pooled") if seg == "distal" else ("proximal", a.chrom)
        g = M.get(key)
        if g is None:
            continue
        n = len(bl)
        for k, (f, s, e) in enumerate(bl):
            if f not in alphabet:
                continue
            m = list(sy); m[k] = AS.WILDCARD
            _, asg = g.viterbi(m)
            pr, po, ms = predict_slot(g, asg[k], alphabet=alphabet)
            nb = [bl[t][0] for t in range(k - 1, -1, -1) if bl[t][0] in alphabet]
            rows.append(dict(arm=a.id, chrom=a.chrom, source=a.source, seg=seg,
                             truth=f, bp=e - s, pred=pr, post=po,
                             base_major=maj,
                             base_pos=posmod.get((seg, min(int(10 * k / max(n, 1)), 9)), maj),
                             base_neigh=(nb[0] if nb else maj)))
    return pd.DataFrame(rows)


def knn(D, lab, trm, tem, k=5):
    Dt = D[np.ix_(tem, trm)]; ytr = lab[trm]; pred = []
    for r in range(Dt.shape[0]):
        idx = np.argsort(Dt[r])[:k]
        w = 1.0 / (Dt[r][idx] + 1e-6)
        sc = collections.Counter()
        for i, ww in zip(idx, w):
            sc[ytr[i]] += ww
        pred.append(sc.most_common(1)[0][0])
    return np.array(pred), lab[tem]


# ---------------------------------------------------------------- main
def main(work, out, fig):
    os.makedirs(out, exist_ok=True); os.makedirs(fig, exist_ok=True)
    rng = np.random.default_rng(SEED)

    ped = AS.load_pedigree(os.path.join(work, "all_samples_censat.bed"))
    hp = AS.load_censat_dir(os.path.join(work, "hprc"), os.path.join(work, "aux"))
    arms = ped + hp
    print(f"[1] arms: {len(arms)}  ({len(ped)} pedigree, {len(hp)} HPRC from "
          f"{len({a.sample for a in hp})} individuals)")

    pd.DataFrame([dict(id=a.id, sample=a.sample, hap=a.hap, chrom=a.chrom,
                       source=a.source, arm_mb=a.length / 1e6, n_blocks=len(a.blocks),
                       n_unspecified=sum(1 for b in a.blocks if b[0] == AS.WILDCARD),
                       unspecified_mb=sum(b[2] - b[1] for b in a.blocks
                                          if b[0] == AS.WILDCARD) / 1e6,
                       syntax="-".join(a.collapsed())) for a in arms]
                 ).to_csv(os.path.join(out, "arm_inventory.csv"), index=False)

    segs = {}
    for a in arms:
        s = AS.split_segments(a)
        if s:
            segs[a.id] = (a, s[0], s[1])
    print(f"[2] arms with an rDNA anchor: {len(segs)}")

    # -- grammars fitted on everything, for description and for placement ----
    M_all, maj_all, pos_all = fit_models(segs, list(segs), rounds=12)
    rows = []
    for (seg, c), g in M_all.items():
        for j in range(g.L):
            top = np.argsort(-g.E[j])[:3]
            rows.append(dict(group=c, segment=seg, slot=j, seed=g.slots[j],
                             **{f"top{k+1}": f"{g.syms[t]}:{g.E[j,t]:.3f}"
                                for k, t in enumerate(top)}))
    pd.DataFrame(rows).to_csv(os.path.join(out, "slot_grammar_emissions.csv"), index=False)

    # -- benchmark: held-out individuals, and cross-dataset ------------------
    samples = sorted({segs[i][0].sample for i in segs})
    rng.shuffle(samples); cut = int(0.7 * len(samples))
    tr_s = set(samples[:cut])
    tr_ids = [i for i in segs if segs[i][0].sample in tr_s]
    te_ids = [i for i in segs if segs[i][0].sample not in tr_s]
    ped_ids = [i for i in segs if segs[i][0].source == "pedigree"]
    hpr_ids = [i for i in segs if segs[i][0].source == "censat"]

    bench = []
    for name, tri, tei in [("held-out individuals (70/30)", tr_ids, te_ids),
                           ("train pedigree -> test HPRC", ped_ids, hpr_ids),
                           ("train HPRC -> test pedigree", hpr_ids, ped_ids)]:
        M, mj, ps = fit_models(segs, tri)
        E = mask_and_recover(segs, M, mj, ps, tei, rng)
        if not len(E):
            continue
        nm = E[E.pred.notna()]
        bench.append(dict(set=name, blocks=len(E),
                          named_pct=100 * E.pred.notna().mean(),
                          acc_among_named=100 * (nm.pred == nm.truth).mean() if len(nm) else np.nan,
                          acc_overall=100 * (E.pred.fillna("_") == E.truth).mean(),
                          base_position=100 * (E.base_pos == E.truth).mean(),
                          base_majority=100 * (E.base_major == E.truth).mean(),
                          base_neighbour=100 * (E.base_neigh == E.truth).mean()))
        if name.startswith("held-out"):
            E.to_csv(os.path.join(out, "placement_eval_heldout.csv"), index=False)
            cal = E[E.pred.notna()].copy()
            cal["bin"] = pd.cut(cal.post, [0, .5, .7, .8, .9, .95, .99, 1.001], right=False)
            (cal.assign(correct=cal.pred == cal.truth)
                .groupby("bin", observed=True)
                .agg(n=("correct", "size"), acc=("correct", "mean"))
                .to_csv(os.path.join(out, "placement_calibration.csv")))
    B = pd.DataFrame(bench)
    B.to_csv(os.path.join(out, "placement_benchmark.csv"), index=False)
    print("[3] masked-block recovery\n", B.round(1).to_string(index=False))

    # -- landmark recovery: labels the CenSat vocabulary does not carry ------
    rows = []
    for i in ped_ids:
        a, d, p = segs[i]
        for seg, bl in (("distal", d), ("proximal", p)):
            key = ("distal", "pooled") if seg == "distal" else ("proximal", a.chrom)
            g = M_all.get(key)
            if g is None or len(bl) < 4:
                continue
            sy = [b[0] for b in bl]
            for k, (f, s, e) in enumerate(bl):
                if f not in ("ACRO", "SST1", "DJ", "PJ"):
                    continue
                m = list(sy); m[k] = AS.WILDCARD
                _, asg = g.viterbi(m)
                pr, po, _ = predict_slot(g, asg[k], alphabet=SATS2 + ["DJ", "PJ"])
                rows.append(dict(arm=a.id, chrom=a.chrom, seg=seg, truth=f,
                                 pred=pr, post=po, bp=e - s))
    LV = pd.DataFrame(rows)
    LV.to_csv(os.path.join(out, "landmark_recovery.csv"), index=False)

    # -- place the real unspecified-satellite blocks -------------------------
    up = []
    for i, (a, d, p) in segs.items():
        for seg, bl in (("distal", d), ("proximal", p)):
            key = ("distal", "pooled") if seg == "distal" else ("proximal", a.chrom)
            g = M_all.get(key)
            if g is None or len(bl) < 4:
                continue
            sy = [b[0] for b in bl]
            _, asg = g.viterbi(sy)
            for (f, s, e), j in zip(bl, asg):
                if f != AS.WILDCARD:
                    continue
                pr, po, ms = predict_slot(g, j)
                up.append(dict(arm=a.id, chrom=a.chrom, source=a.source, segment=seg,
                               slot=j, start=s, end=e, bp=e - s, pred=pr, post=po,
                               slot_sat_mass=ms))
    U = pd.DataFrame(up)
    U.to_csv(os.path.join(out, "unknown_block_placements.csv"), index=False)
    tot = {c: sum(1 for v in segs.values() if v[0].chrom == c) for c in CHR}
    G = (U[U.slot.ge(0)].groupby(["chrom", "segment", "slot"])
         .agg(arms=("arm", "nunique"), Mb=("bp", lambda s: s.sum() / 1e6),
              med_kb=("bp", lambda s: s.median() / 1e3),
              named_frac=("pred", lambda s: s.notna().mean()),
              top_pred=("pred", lambda s: s.value_counts().index[0] if s.notna().any() else None),
              med_post=("post", "median")).reset_index())
    G["pct_arms"] = [100 * r.arms / tot[r.chrom] for r in G.itertuples()]
    G.to_csv(os.path.join(out, "unknown_slot_summary.csv"), index=False)
    print(f"[4] unspecified-satellite blocks: {len(U)}; "
          f"{U.slot.ge(0).sum()} placed, {U.pred.notna().sum()} named")

    # -- chromosome assignment from syntax alone ----------------------------
    lab = np.array([a.chrom for a in arms]); src = np.array([a.source for a in arms])
    smp = np.array([a.sample for a in arms])
    us = sorted(set(smp)); rng.shuffle(us)
    trm = np.isin(smp, us[:int(0.7 * len(us))]); tem = ~trm
    ctl = []
    for wt, orders, name in [("sqrt_bp", (1, 2, 3), "size-weighted n-grams"),
                             ("none", (1, 2, 3), "order only, unweighted"),
                             ("none", (2, 3), "order only, 2+3-grams"),
                             ("sqrt_bp", (1,), "composition only")]:
        X, _ = AS.ngram_features(arms, orders=orders, weight=wt)
        D = AS.cosine_distance(X)
        p, y = knn(D, lab, trm, tem)
        ctl.append(dict(representation=name, accuracy_pct=100 * (p == y).mean()))
        if name == "size-weighted n-grams":
            pd.crosstab(pd.Series(y, name="truth"), pd.Series(p, name="assigned")
                        ).to_csv(os.path.join(out, "chrom_assignment_confusion.csv"))
            for nm2, a_, b_ in [("train pedigree -> test HPRC", src == "pedigree", src == "censat"),
                                ("train HPRC -> test pedigree", src == "censat", src == "pedigree")]:
                p2, y2 = knn(D, lab, a_, b_)
                ctl.append(dict(representation=nm2, accuracy_pct=100 * (p2 == y2).mean()))
    Ln = np.array([[a.length] for a in arms], float)
    p, y = knn(np.abs(Ln - Ln.T) / 1e6, lab, trm, tem)
    ctl.append(dict(representation="arm length only", accuracy_pct=100 * (p == y).mean()))
    ctl.append(dict(representation="chance (largest class)",
                    accuracy_pct=100 * pd.Series(lab).value_counts(normalize=True).max()))
    pd.DataFrame(ctl).to_csv(os.path.join(out, "chrom_assignment_controls.csv"), index=False)
    print("[5] chromosome assignment\n", pd.DataFrame(ctl).round(1).to_string(index=False))
    print("\ndone; tables in", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="../work")
    ap.add_argument("--out", default="../data")
    ap.add_argument("--fig", default="../figures")
    a = ap.parse_args()
    main(a.work, a.out, a.fig)
