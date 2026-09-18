"""
repeatcn - copy number and array mass of multi-copy sequence classes from WGS.

Scope. This estimates how much of a repeat class an individual carries. It does not
assemble, phase, or locate it. For the acrocentric short arms that is often the only
measurement short reads can support (see report sections 19 and 22.5), and for rDNA it
is a measurement with demonstrated phenotypic associations.

Design. Depth extraction is delegated to mosdepth or samtools, which are fast and
correct; this module owns the parts that are easy to get wrong:

  * the estimator and its ploidy/read-length bookkeeping
  * choice of denominator, including the split denominator needed whenever two classes
    are compared with each other (see `estimate(..., denominator='split')`)
  * GC-bias correction fitted from the sample's own single-copy windows
  * cohort technical covariates taken from coverage PCs rather than from batch labels,
    following NGS-PCA (https://github.com/jlanej/NGS-PCA)
  * negative-control classes and within-class concordance as bias detectors

Estimator. Write d for depth against the haploid reference, so a "30x genome" has
d = 30 and a diploid single-copy locus also reads d = 30. For a class c with observed
base count B_c:

    diploid array mass   M_c  = 2 * B_c / d_ctrl                      [bp]
    diploid copy number  CN_c = 2 * B_c / (d_ctrl * L_unit)           [copies]

where d_ctrl is the GC-corrected depth over single-copy autosomal control windows.
Counting in k-mer space instead of base space multiplies numerator and denominator by
the same (R - k + 1)/R factor, so the same formulae apply unchanged to k-mer counts -
provided the control is counted the same way.

Which quantity to report is a property of the class, not a preference. Classes with a
defined repeating unit (rDNA, 5S, D4Z4, DXZ4, SST1) get copy number. Classes whose
higher-order structure varies between arrays and individuals (HSat1, HSat2, HSat3) have
no stable unit, and reporting "copies" for them invents precision that does not exist;
they get array mass in Mb.

Subcommands:
    targets    build the target BED and class table from CHM13 annotation
    estimate   turn mosdepth/samtools region depths into CN and mass estimates
    adjust     residualise a cohort's estimates on NGS-PCA coverage PCs
    selftest   validate the estimator, the artefacts it avoids and the PC adjustment

Usage:
    python repeatcn.py targets  --annot ../data --out ../data
    mosdepth --by targets.bed --no-per-base --fast-mode sample sample.cram
    python repeatcn.py estimate --regions sample.regions.bed.gz --targets targets.tsv \\
                                --out sample.cn.tsv
    python repeatcn.py selftest
"""
import argparse
import csv
import gzip
import io
import math
import os
import random
import statistics
import sys
from collections import defaultdict

# ----------------------------------------------------------------- class table

# Unit lengths. "computed" values are measured from a reference record in this repo;
# "literature" values are the conventional consensus unit sizes for the class. A class
# with unit=None has no stable unit and is reported as array mass.
UNITS = {
    "rDNA45S":  (44838, "computed",   "KY962518.1 complete repeating unit"),
    "rDNA5S":   (2213,  "literature", "1q42 tandem array unit"),
    "D4Z4":     (3300,  "literature", "4q35 and 10q26 macrosatellite; FSHD locus"),
    "DXZ4":     (3000,  "literature", "Xq23 macrosatellite"),
    "SST1":     (1400,  "literature", "NBL2 macrosatellite"),
    "alphaSat": (171,   "literature", "alpha-satellite monomer"),
    "betaSat":  (68,    "literature", "Sau3A/BSR monomer"),
    "gammaSat": (220,   "literature", "gamma-satellite monomer"),
    "mtDNA":    (16569, "computed",   "NC_012920.1"),
    "TELO":     (6,     "literature", "TTAGGG; prefer the k-mer backend, reference "
                                      "telomeres are not representative of sample length"),
    "HSat1A":   (None,  "-",          "no stable higher-order unit"),
    "HSat1B":   (None,  "-",          "no stable higher-order unit"),
    "HSat2":    (None,  "-",          "pentamer-derived, variable HOR"),
    "HSat3":    (None,  "-",          "pentamer-derived, variable HOR"),
    "ACRO":     (None,  "-",          "acrocentric composite repeat"),
    "LSAU-BSAT": (None, "-",          "composite"),
}

# CenSat v2.1 name prefix -> class. 'ct' (non-satellite transition) and the generic
# 'censat' bucket are deliberately excluded: they are not a repeat class.
CENSAT_MAP = {"hsat1A": "HSat1A", "hsat1B": "HSat1B", "hsat2": "HSat2", "hsat3": "HSat3",
              "bsat": "betaSat", "gsat": "gammaSat", "hor": "alphaSat_HOR",
              "dhor": "alphaSat_dHOR", "mon": "alphaSat_mon", "rDNA": "rDNA45S"}
COMPOSITE_MAP = {"5SrRNA": "rDNA5S", "ACRO": "ACRO", "LSAU-BSAT": "LSAU-BSAT",
                 "TELO": "TELO"}

# Classes whose reference annotation is absent from the CHM13 satellite tracks used
# here; they need a panel built from their own consensus before they can be estimated.
NEEDS_PANEL = ["D4Z4", "DXZ4", "mtDNA"]


def read_bed(path, name_col=3):
    rows = []
    op = gzip.open if path.endswith(".gz") else io.open
    with op(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > name_col and p[1].isdigit():
                rows.append((p[0], int(p[1]), int(p[2]), p[name_col]))
    return rows


def build_targets(annot_dir, outdir, coverage=30.0, read_len=150.0):
    """Assemble the target intervals and the per-class summary table."""
    inter = defaultdict(list)
    cen = os.path.join(annot_dir, "chm13v2.0_censat_v2.1.bed")
    for chrom, s, e, name in read_bed(cen):
        pref = name.split("_")[0].split("(")[0]
        cls = CENSAT_MAP.get(pref)
        if cls:
            inter[cls].append((chrom, s, e))
        if "SST1" in name.upper():
            inter["SST1"].append((chrom, s, e))
    comp = os.path.join(annot_dir, "chm13v2.0_composite-repeats_2022DEC.bed")
    if os.path.exists(comp):
        for chrom, s, e, name in read_bed(comp):
            cls = COMPOSITE_MAP.get(name.rsplit("_", 1)[0])
            if cls:
                inter[cls].append((chrom, s, e))

    rows = []
    for cls, ivs in sorted(inter.items()):
        bp = sum(e - s for _, s, e in ivs)
        base = cls.split("_")[0]
        unit, prov, note = UNITS.get(base, (None, "-", ""))
        reads = coverage * bp / read_len          # haploid mass, so a per-haplotype count
        pct_se = 100.0 / math.sqrt(reads) if reads else float("nan")
        rows.append(dict(
            cls=cls, quantity="copies" if unit else "megabases",
            unit_bp=unit or "", unit_source=prov,
            chm13_haploid_Mb=round(bp / 1e6, 3),
            implied_haploid_copies=round(bp / unit) if unit else "",
            n_arrays=len(ivs), n_chrom=len({c for c, _, _ in ivs}),
            median_array_kb=round(statistics.median([e - s for _, s, e in ivs]) / 1e3, 1),
            reads_at_30x=int(reads), sampling_pct_SE=round(pct_se, 3), note=note))
    for cls in NEEDS_PANEL:
        unit, prov, note = UNITS[cls]
        rows.append(dict(cls=cls, quantity="copies", unit_bp=unit, unit_source=prov,
                         chm13_haploid_Mb="", implied_haploid_copies="", n_arrays="",
                         n_chrom="", median_array_kb="", reads_at_30x="",
                         sampling_pct_SE="", note=note + "; not in the CHM13 satellite "
                         "tracks used here - supply a consensus panel"))
    rows.sort(key=lambda r: -(r["chm13_haploid_Mb"] or 0))

    os.makedirs(outdir, exist_ok=True)
    tsv = os.path.join(outdir, "repeat_cn_targets.tsv")
    with io.open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                            lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    bed = os.path.join(outdir, "repeat_cn_targets.bed")
    with io.open(bed, "w") as fh:
        for cls, ivs in sorted(inter.items()):
            for chrom, s, e in sorted(ivs):
                fh.write("%s\t%d\t%d\t%s\n" % (chrom, s, e, cls))
    return rows, tsv, bed


# ----------------------------------------------------------------- GC correction

def gc_curve(windows, nbins=40, min_windows=50):
    """Median depth by GC bin over single-copy control windows, normalised to 1.

    windows: iterable of (gc_fraction, depth). Returns a function f(gc) -> factor.
    Bins with too few windows fall back to 1.0 rather than to a noisy estimate.
    """
    by = defaultdict(list)
    for gc, d in windows:
        by[min(nbins - 1, max(0, int(gc * nbins)))].append(d)
    alld = [d for _, d in windows]
    if not alld:
        return lambda gc: 1.0
    overall = statistics.median(alld)
    if not overall:
        return lambda gc: 1.0
    # one knot per sufficiently populated bin, placed at the bin's median GC rather
    # than its centre, so a bin that is unevenly filled is not mis-located
    knots = sorted((statistics.median([g for g, _ in windows
                                       if min(nbins - 1, int(g * nbins)) == b]),
                    statistics.median(v) / overall)
                   for b, v in by.items() if len(v) >= min_windows)
    if not knots:
        return lambda gc: 1.0

    def f(gc):
        # linear interpolation between knots; flat extrapolation outside their range.
        # Satellite classes sit far from the genomic mean GC, where the bias curve is
        # steepest, so nearest-bin lookup biases the correction by several percent.
        if gc <= knots[0][0]:
            return knots[0][1]
        if gc >= knots[-1][0]:
            return knots[-1][1]
        for i in range(1, len(knots)):
            g1, f1 = knots[i - 1]
            g2, f2 = knots[i]
            if gc <= g2:
                w = (gc - g1) / (g2 - g1) if g2 > g1 else 0.0
                return f1 + w * (f2 - f1)
        return knots[-1][1]
    return f


# ----------------------------------------------------------------- estimator

def estimate(class_bases, class_ref_bp, control_depth, unit_bp=None,
             class_gc=None, gc_f=None):
    """Diploid array mass (bp) and, where a unit is defined, diploid copy number.

    class_bases   sequenced bases attributed to the class
    class_ref_bp  reference span the depth was measured over (unused for mass, kept
                  for callers that want the implied depth back)
    control_depth GC-corrected depth over single-copy autosomal control windows,
                  expressed against the haploid reference
    """
    if control_depth <= 0:
        raise ValueError("control depth must be positive")
    d = control_depth
    if gc_f is not None and class_gc is not None:
        d = d * gc_f(class_gc)          # expected depth at the class's GC
    mass = 2.0 * class_bases / d
    out = dict(diploid_mass_bp=mass, diploid_mass_Mb=mass / 1e6,
               effective_control_depth=d)
    if unit_bp:
        out["diploid_copies"] = mass / unit_bp
    if class_ref_bp:
        out["class_depth"] = class_bases / class_ref_bp
    return out


def control_depth_from_windows(windows, denominator="all", seed=0):
    """Single-copy control depth. 'split' returns two independent halves.

    Any two classes normalised by the *same* denominator inherit a positive correlation
    from the denominator's own noise: in logs the shared term contributes variance
    sigma_n^2 to the covariance. When the analysis correlates two classes with each
    other - rDNA against mtDNA being the obvious case - normalise them by disjoint
    halves of the control windows so that term is zero.
    """
    d = [x[1] for x in windows]
    if denominator == "all":
        return statistics.median(d)
    if denominator == "split":
        r = random.Random(seed)
        idx = list(range(len(d)))
        r.shuffle(idx)
        h = len(idx) // 2
        return (statistics.median([d[i] for i in idx[:h]]),
                statistics.median([d[i] for i in idx[h:]]))
    raise ValueError("denominator must be 'all' or 'split'")


# ------------------------------------------------- batch, estimated from the data

def coverage_pcs(bin_matrix, n_pc=20):
    """PCs of median-centered log2 fold-change depth, following NGS-PCA.

    bin_matrix: samples x bins array of mean depth per fixed-width bin, from mosdepth
    --by, restricted to autosomal bins that do not overlap the exclusion set (SV
    blacklist, low mappability, DGV, segmental duplications).

    Normalisation follows jlanej/NGS-PCA: within each sample take log2 of depth over
    that sample's median bin depth, then centre each bin to a median of zero across
    samples. Median rather than mean centering makes this a robust variant of PCA, so
    methods text should say "PCA on median-centered log2 fold change of depth". NGS-PCA
    itself uses a randomized SVD (Halko et al. 2011) to make this tractable at cohort
    scale; the exact SVD here is for modest N and for the self-test. For real cohorts
    run NGS-PCA and read its svd.pcs.txt with `load_pcs`.

    Returns (scores samples x n_pc, variance_explained per component).
    """
    import numpy as np
    X = np.asarray(bin_matrix, float)
    med = np.median(X, axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log2(X / med)
    L[~np.isfinite(L)] = 0.0
    L = L - np.median(L, axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(L, full_matrices=False)
    k = min(n_pc, len(S))
    return U[:, :k] * S[:k], (S ** 2 / (S ** 2).sum())[:k]


def load_auto_hq_median(path, strip_suffix=None):
    """Read NGS-PCA's autosomal.median.txt -> {sample: (median_depth, n_bins)}.

    Columns are SAMPLE, AUTO_HQ_median, N_BINS. This is the preferred denominator for
    the estimator: it is the per-sample median depth over exactly the bins the PCA
    used - autosomal, with the exclusion set (SV blacklist, low mappability, DGV,
    segmental duplications) removed - so it is a single-copy denominator by
    construction, and it is already computed, with no extra pass over the data.

    Two things to know. Values are not floored, so a failed or empty sample reads as
    zero; those are dropped here and returned in the second element. And sample names
    are the mosdepth file name minus 'regions.bed.gz', which leaves a trailing dot
    unless the run used -sampleSuffix; pass strip_suffix to match your own ids.

    Returns (kept {sample: (median, n_bins)}, failed [sample, ...]).
    """
    kept, failed = {}, []
    with io.open(path) as fh:
        head = fh.readline()
        if "AUTO_HQ_median" not in head:
            raise ValueError("%s does not look like autosomal.median.txt" % path)
        for line in fh:
            p = line.split()
            if len(p) < 3:
                continue
            name, med, nb = p[0], float(p[1]), int(p[2])
            if strip_suffix and name.endswith(strip_suffix):
                name = name[: -len(strip_suffix)]
            (failed.append(name) if med <= 0 else
             kept.__setitem__(name, (med, nb)))
    return kept, failed


def median_sampling_se(n_bins, bin_cv=0.20):
    """Relative standard error of the high-quality median, normal approximation.

    Used to check that the denominator's own sampling noise is small enough to ignore
    when two classes share it (see the induced-correlation formula in the module
    docstring for `control_depth_from_windows`). At 15,000 retained bins this is about
    0.2%, which induces a correlation of 5e-5 between two classes - four orders of
    magnitude below what matters. The residual risk from a shared denominator is
    therefore systematic, not statistical, which is what the coverage PCs address.
    """
    return 1.2533 * bin_cv / math.sqrt(max(1, n_bins))


def load_pcs(path, n_pc=None):
    """Read an NGS-PCA svd.pcs.txt (samples in rows, PCs in columns)."""
    pcs = {}
    with io.open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ncol = len(head) - 1
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                pcs[p[0]] = [float(x) for x in p[1:1 + (n_pc or ncol)]]
    return pcs


def adjust_for_pcs(values, pcs, n_pc=None, log=True):
    """Residualise per-sample estimates on coverage PCs.

    values: {sample: estimate}. pcs: {sample: [pc1, ...]}.

    The PC basis is built from single-copy bins with satellite, segmental duplications
    and low-mappability sequence excluded, so it is computed on sequence disjoint from
    every class in the target table. That disjointness is what makes this safe: the
    components can absorb library, batch and GC-bias structure, but they cannot absorb
    the dosage of the class being measured. The residual is returned on the original
    scale with the cohort mean restored.

    Returns (adjusted {sample: value}, variance removed, n_pc used).
    """
    import numpy as np
    ids = sorted(set(values) & set(pcs))
    if len(ids) < 10:
        raise ValueError("need at least 10 samples with both an estimate and PCs")
    y = np.array([values[s] for s in ids], float)
    if log:
        if (y <= 0).any():
            raise ValueError("log adjustment needs positive estimates")
        y = np.log(y)
    P = np.array([pcs[s][:n_pc] if n_pc else pcs[s] for s in ids], float)
    X = np.column_stack([np.ones(len(ids)), P])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    sst = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / sst if sst else 0.0
    out = resid + y.mean()
    return dict(zip(ids, np.exp(out) if log else out)), r2, X.shape[1] - 1


# ------------------------------------------------- validation against trios

def trio_reliability(trios, values, log=True):
    """Estimate the reliability of a CN measurement from parent-offspring trios.

    Array copy number is a physical DNA quantity, so it is inherited additively with
    no dominance and no environment: a child's diploid CN is one paternal haplotype
    plus one maternal haplotype, and therefore

        E[T_child | parents] = (T_father + T_mother) / 2

    holds *exactly*, with a coefficient of exactly 1 - unlike a phenotype, where the
    midparent slope estimates heritability. Regressing the observed child value on the
    observed midparent value therefore has expected slope

        slope = var_true / (var_true + var_err) = reliability

    so the departure of the slope from 1 is measurement error, and nothing else. With
    the observed parental variance V this gives the full decomposition directly:
    var_true = R * V and var_err = (1 - R) * V.

    Mutation does not bias the slope - it is mean-zero and independent of the parents,
    so it adds variance to the child only. It appears as residual variance in excess
    of the model's prediction, which is reported as `excess_resid_cv`.

    **The one thing that breaks this:** measurement error shared within a family. If a
    trio was extracted, libraried and sequenced together, part of e is common to all
    three, which enters the numerator and pushes the slope up. With all of the error
    shared the slope is exactly 1 however bad the assay is. Run this before and after
    `adjust_for_pcs`; if reliability falls after adjustment, the raw figure was
    inflated by shared batch. Never report the unadjusted number alone.

    trios: [(father_id, mother_id, child_id)]. values: {sample: estimate}.
    """
    import numpy as np
    f, m, c = [], [], []
    for fa, mo, ch in trios:
        if fa in values and mo in values and ch in values:
            v = [values[fa], values[mo], values[ch]]
            if log and min(v) <= 0:
                continue
            v = [math.log(x) for x in v] if log else v
            f.append(v[0]); m.append(v[1]); c.append(v[2])
    n = len(c)
    if n < 30:
        raise ValueError("need at least 30 complete trios, got %d" % n)
    f, m, c = np.array(f), np.array(m), np.array(c)
    mp = (f + m) / 2.0
    b, a = np.polyfit(mp, c, 1)
    resid = c - (a + b * mp)
    se_b = float(np.sqrt((resid ** 2).sum() / (n - 2) / ((mp - mp.mean()) ** 2).sum()))
    V = float(np.var(np.concatenate([f, m]), ddof=1))       # observed parental variance
    R = float(min(max(b, 0.0), 1.0))
    var_true, var_err = R * V, (1 - R) * V
    pred_resid = V * (1 - R ** 2 / 2)
    obs_resid = float(np.var(resid, ddof=2))
    excess = max(0.0, obs_resid - pred_resid)
    return dict(n_trios=n, reliability=round(float(b), 4), se=round(se_b, 4),
                reliability_ci=(round(float(b) - 1.96 * se_b, 4),
                                round(float(b) + 1.96 * se_b, 4)),
                total_cv=round(100 * math.sqrt(V), 2),
                true_cv=round(100 * math.sqrt(var_true), 2),
                error_cv=round(100 * math.sqrt(max(var_err, 0)), 2),
                resid_observed=round(obs_resid, 5), resid_predicted=round(pred_resid, 5),
                excess_resid_cv=round(100 * math.sqrt(excess), 2))


def load_trios(path):
    """Read trios from a PLINK-style .ped/.fam (IID, PAT, MAT in columns 2-4) or a
    three-column father/mother/child file. Returns [(father, mother, child)]."""
    trios = []
    with io.open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 6 and p[2] not in ("0", "-9", "") and p[3] not in ("0", "-9", ""):
                trios.append((p[2], p[3], p[1]))          # PAT, MAT, IID
            elif len(p) == 3 and p[0].lower() not in ("father", "fa", "paternal"):
                trios.append((p[0], p[1], p[2]))
    return trios


def permuted_trio_null(trios, values, n_perm=200, seed=0, log=True):
    """Reliability after shuffling children between families: should be ~0.

    A non-zero null means the estimate carries structure shared across unrelated
    samples - batch, ancestry, or depth - rather than transmitted biology.
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    kids = [t[2] for t in trios]
    out = []
    for _ in range(n_perm):
        perm = list(kids)
        rng.shuffle(perm)
        try:
            out.append(trio_reliability([(t[0], t[1], k) for t, k in zip(trios, perm)],
                                        values, log)["reliability"])
        except ValueError:
            pass
    return dict(mean=round(float(np.mean(out)), 4), sd=round(float(np.std(out)), 4),
                n_perm=len(out))


def concordance(subregion_estimates):
    """Coefficient of variation across sub-regions of one class.

    For rDNA, estimate separately over 18S, 5.8S, 28S and the IGS. The subunits are
    physically the same molecule at the same copy number, so disagreement between them
    measures bias, not biology. Report it; do not average it away.
    """
    v = [x for x in subregion_estimates if x and x > 0]
    if len(v) < 2:
        return None
    m = statistics.mean(v)
    return statistics.pstdev(v) / m if m else None


def estimate_from_regions(regions_path, targets_path, control_gc=None,
                          denominator="all", control_depth=None):
    """Read mosdepth --by output and a targets table; return per-class estimates.

    mosdepth writes <prefix>.regions.bed.gz with columns chrom, start, end, name, mean
    depth. Control windows are supplied separately (a mosdepth run over single-copy
    autosomal windows with a GC column), because the right control set depends on the
    reference and the callset.
    """
    tgt = {}
    with io.open(targets_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tgt[r["cls"]] = r
    bases, span = defaultdict(float), defaultdict(float)
    op = gzip.open if regions_path.endswith(".gz") else io.open
    with op(regions_path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            s, e, name, dep = int(p[1]), int(p[2]), p[3], float(p[4])
            bases[name] += dep * (e - s)
            span[name] += (e - s)
    if control_depth is None and control_gc is None:
        raise ValueError("supply control_depth (NGS-PCA AUTO_HQ_median, preferred) "
                         "or control windows; see docs/REPEAT_CN.md")
    # GC correction still needs windows; the scale can come from the median alone
    gcf = gc_curve(control_gc) if control_gc else (lambda gc: 1.0)
    ctrl = (control_depth if control_depth is not None
            else control_depth_from_windows(control_gc, denominator))
    out = []
    for cls, b in sorted(bases.items()):
        r = tgt.get(cls, {})
        unit = int(r["unit_bp"]) if r.get("unit_bp") else None
        c = ctrl[0] if isinstance(ctrl, tuple) else ctrl
        est = estimate(b, span[cls], c, unit_bp=unit, gc_f=gcf)
        est.update(cls=cls, quantity=r.get("quantity", ""),
                   reference_Mb=r.get("chm13_haploid_Mb", ""))
        out.append(est)
    return out


# ----------------------------------------------------------------- self-test

def selftest(n=200000, seed=7):
    """Validate the estimator and the two failure modes it is built to avoid.

    No data required: counts are simulated with known truth.
    """
    rng = random.Random(seed)
    ok = True

    # 1. recovery of a known copy number
    true_cn, unit, depth = 380.0, 44838, 30.0
    bases = true_cn / 2 * unit * depth          # haploid copies x unit x depth
    got = estimate(bases, true_cn / 2 * unit, depth, unit_bp=unit)["diploid_copies"]
    ok &= abs(got - true_cn) < 1e-6
    print("[1] copy-number recovery      : %.3f vs %.3f  %s"
          % (got, true_cn, "OK" if abs(got - true_cn) < 1e-6 else "FAIL"))

    # 2. shared vs split denominator, correlating two independent classes
    sr, sm, sn = 0.25, 0.30, 0.05               # log SDs: class1, class2, denominator
    shared, split = [], []
    for _ in range(n):
        er, em, en = rng.gauss(0, sr), rng.gauss(0, sm), rng.gauss(0, sn)
        en2 = rng.gauss(0, sn)
        shared.append((er - en, em - en))
        split.append((er - en, em - en2))

    def pearson(p):
        xs, ys = [a for a, _ in p], [b for _, b in p]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        cov = sum((x - mx) * (y - my) for x, y in p)
        return cov / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    r_shared, r_split = pearson(shared), pearson(split)
    predicted = sn ** 2 / math.sqrt((sr ** 2 + sn ** 2) * (sm ** 2 + sn ** 2))
    tol = 4.0 / math.sqrt(n)                    # 4 standard errors of a correlation
    p2 = abs(r_shared - predicted) < tol and abs(r_split) < tol
    ok &= p2
    print("[2] denominator artefact      : shared r = %+.4f (predicted %+.4f), "
          "split r = %+.4f, tol %.4f  %s"
          % (r_shared, predicted, r_split, tol, "OK" if p2 else "FAIL"))

    # 3. GC correction: a GC-rich class under a simulated bias curve
    def bias(gc):                               # positive everywhere, peaks at 0.41
        return math.exp(-((gc - 0.41) / 0.13) ** 2)
    wins = [((g := rng.uniform(0.30, 0.65)), 30.0 * bias(g) * rng.gauss(1, 0.02))
            for _ in range(20000)]
    f = gc_curve(wins)
    ctrl = statistics.median([d for _, d in wins])
    class_gc, true_mass = 0.581, 17.0e6         # rDNA-like: 58.1% GC, 17 Mb diploid
    observed = true_mass / 2 * 30.0 * bias(class_gc)
    e_naive = 100 * (estimate(observed, true_mass / 2, ctrl)["diploid_mass_bp"]
                     - true_mass) / true_mass
    e_corr = 100 * (estimate(observed, true_mass / 2, ctrl, class_gc=class_gc,
                             gc_f=f)["diploid_mass_bp"] - true_mass) / true_mass
    p3 = abs(e_corr) < 2.0
    ok &= p3
    print("[3] GC correction (58%% GC)    : error %+.1f%% -> %+.1f%%  %s"
          % (e_naive, e_corr, "OK" if p3 else "FAIL"))

    # 4. concordance flags a biased sub-region
    c_clean, c_bias = concordance([380, 379, 381, 380]), concordance([380, 379, 381, 300])
    p4 = c_clean < 0.01 < c_bias
    ok &= p4
    print("[4] subunit concordance       : clean %.3f, one biased subunit %.3f  %s"
          % (c_clean, c_bias, "OK" if p4 else "FAIL"))
    # 5. coverage PCs remove batch confounding without removing biology
    try:
        import numpy as np
        ns, nb, nbatch = 400, 300, 4
        rs = np.random.default_rng(seed)
        batch = rs.integers(0, nbatch, ns)
        sig = rs.normal(0, 1.0, (nbatch, nb))            # each batch's bin signature
        shift_cls = rs.normal(0, 0.30, nbatch)           # batch effect on the class
        shift_phe = rs.normal(0, 0.80, nbatch)           # and on the phenotype
        true_log = rs.normal(0, 0.25, ns)                # biology
        pheno = (0.15 * true_log / 0.25 + shift_phe[batch] + rs.normal(0, 1.0, ns))
        obs = np.exp(true_log + shift_cls[batch] + rs.normal(0, 0.02, ns))
        bins = 30.0 * np.exp(0.10 * sig[batch] + rs.normal(0, 0.02, (ns, nb)))
        sc, _ = coverage_pcs(bins, n_pc=nbatch)
        ids = ["s%d" % i for i in range(ns)]
        adj, r2, npc = adjust_for_pcs(dict(zip(ids, obs)),
                                      {s: list(sc[i]) for i, s in enumerate(ids)})
        adj_v = np.array([adj[s] for s in ids])

        def cor(a, b):
            a, b = np.asarray(a, float), np.asarray(b, float)
            return float(np.corrcoef(a, b)[0, 1])
        r_true = cor(true_log, pheno)
        r_raw = cor(np.log(obs), pheno)
        r_adj = cor(np.log(adj_v), pheno)
        p5 = abs(r_adj - r_true) < abs(r_raw - r_true) / 2
        ok &= p5
        print("[5] coverage-PC adjustment    : r(pheno) true %+.3f, raw %+.3f, "
              "adjusted %+.3f; %d PCs remove %.0f%% of estimate variance  %s"
              % (r_true, r_raw, r_adj, npc, 100 * r2, "OK" if p5 else "FAIL"))
    except ImportError:
        print("[5] coverage-PC adjustment    : skipped (numpy not available)")

    # 6. trio reliability, with and without family-shared error
    try:
        import numpy as np
        rs = np.random.default_rng(seed + 1)
        nt, s_hap, s_err = 602, 0.18, 0.12
        R_true = 2 * s_hap ** 2 / (2 * s_hap ** 2 + s_err ** 2)
        def cohort(shared_frac):
            vals, trios = {}, []
            for i in range(nt):
                hap = rs.normal(0, s_hap, 4)                # fa1 fa2 mo1 mo2
                T = dict(f=hap[0] + hap[1], m=hap[2] + hap[3],
                         c=hap[rs.integers(0, 2)] + hap[2 + rs.integers(0, 2)])
                fam = rs.normal(0, s_err * math.sqrt(shared_frac))
                ids = {}
                for who in "fmc":
                    sid = "t%d_%s" % (i, who)
                    e = fam + rs.normal(0, s_err * math.sqrt(1 - shared_frac))
                    vals[sid] = math.exp(T[who] + e)
                    ids[who] = sid
                trios.append((ids["f"], ids["m"], ids["c"]))
            return trios, vals
        tr0, v0 = cohort(0.0)
        r0 = trio_reliability(tr0, v0)
        tr1, v1 = cohort(1.0)
        r1 = trio_reliability(tr1, v1)
        null = permuted_trio_null(tr0, v0, n_perm=60, seed=3)
        p6 = (abs(r0["reliability"] - R_true) < 3 * r0["se"]
              and r1["reliability"] > R_true + 2 * r1["se"]
              and abs(null["mean"]) < 0.06)
        ok &= p6
        print("[6] trio reliability          : independent error R = %.3f +/- %.3f "
              "(true %.3f); fully shared error R = %.3f; permuted null %+.3f  %s"
              % (r0["reliability"], r0["se"], R_true, r1["reliability"],
                 null["mean"], "OK" if p6 else "FAIL"))
    except ImportError:
        print("[6] trio reliability          : skipped (numpy not available)")

    print("\n%s" % ("all checks passed" if ok else "FAILURES - see above"))
    return ok


# ----------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("targets", help="build target BED and class table")
    t.add_argument("--annot", default="../data")
    t.add_argument("--out", default="../data")
    e = sub.add_parser("estimate", help="estimate CN from mosdepth regions")
    e.add_argument("--regions", required=True)
    e.add_argument("--targets", required=True)
    e.add_argument("--auto-hq-median", help="NGS-PCA autosomal.median.txt (preferred "
                                            "denominator)")
    e.add_argument("--sample", help="sample id within autosomal.median.txt")
    e.add_argument("--control", help="TSV of gc<TAB>depth per window; needed for the "
                                     "GC correction, and as a fallback denominator")
    e.add_argument("--denominator", default="all", choices=["all", "split"])
    e.add_argument("--out", default="-")
    j = sub.add_parser("adjust", help="residualise cohort estimates on coverage PCs")
    j.add_argument("--estimates", required=True,
                   help="TSV: sample<TAB>value (one class), header required")
    j.add_argument("--pcs", required=True, help="NGS-PCA svd.pcs.txt")
    j.add_argument("--n-pc", type=int, default=None)
    j.add_argument("--no-log", action="store_true")
    j.add_argument("--out", default="-")
    v = sub.add_parser("trios", help="reliability of an estimate from parent-offspring "
                                     "trios")
    v.add_argument("--estimates", required=True, help="TSV: sample<TAB>value, header")
    v.add_argument("--pedigree", required=True, help=".ped/.fam or father/mother/child")
    v.add_argument("--pcs", help="NGS-PCA svd.pcs.txt; also reports the adjusted value")
    v.add_argument("--n-pc", type=int, default=20)
    v.add_argument("--perm", type=int, default=200)
    sub.add_parser("selftest", help="validate the estimator with simulated counts")
    a = ap.parse_args(argv)

    if a.cmd == "targets":
        rows, tsv, bed = build_targets(a.annot, a.out)
        print("%d classes -> %s, %s" % (len(rows), tsv, bed))
        for r in rows[:12]:
            print("  %-16s %-9s %8s Mb  %s"
                  % (r["cls"], r["quantity"], r["chm13_haploid_Mb"],
                     ("%s copies" % r["implied_haploid_copies"]) if r["unit_bp"] else ""))
    elif a.cmd == "estimate":
        ctrl = None
        if a.control:
            ctrl = [(float(x[0]), float(x[1])) for x in
                    (l.split() for l in io.open(a.control)) if len(x) >= 2]
        depth = None
        if a.auto_hq_median:
            med, failed = load_auto_hq_median(a.auto_hq_median)
            if not a.sample:
                ap.error("--auto-hq-median requires --sample")
            if a.sample in failed:
                ap.error("sample %s has AUTO_HQ_median 0 - failed or empty" % a.sample)
            if a.sample not in med:
                ap.error("sample %s not in %s (%d present; check -sampleSuffix "
                         "naming)" % (a.sample, a.auto_hq_median, len(med)))
            depth, nbins = med[a.sample]
            sys.stderr.write("denominator: AUTO_HQ_median %.3f over %d bins "
                             "(sampling SE %.3f%%)\n"
                             % (depth, nbins, 100 * median_sampling_se(nbins)))
        res = estimate_from_regions(a.regions, a.targets, ctrl, a.denominator, depth)
        fh = sys.stdout if a.out == "-" else io.open(a.out, "w", newline="")
        # classes without a defined unit carry no diploid_copies field, so the
        # header is the union of keys, not the first row's
        order, seen = [], set()
        for r in res:
            for k in r:
                if k not in seen:
                    seen.add(k)
                    order.append(k)
        for r in res:
            for k in order:
                r.setdefault(k, "")
        w = csv.DictWriter(fh, fieldnames=order, delimiter="\t",
                            lineterminator="\n")
        w.writeheader()
        w.writerows(res)
    elif a.cmd == "adjust":
        vals = {}
        with io.open(a.estimates) as fh:
            fh.readline()
            for line in fh:
                p = line.split()
                if len(p) >= 2:
                    vals[p[0]] = float(p[1])
        adj, r2, npc = adjust_for_pcs(vals, load_pcs(a.pcs, a.n_pc),
                                      n_pc=a.n_pc, log=not a.no_log)
        fh = sys.stdout if a.out == "-" else io.open(a.out, "w")
        fh.write("# %d PCs removed %.1f%% of estimate variance\n" % (npc, 100 * r2))
        fh.write("sample\traw\tadjusted\n")
        for s in sorted(adj):
            fh.write("%s\t%g\t%g\n" % (s, vals[s], adj[s]))
    elif a.cmd == "trios":
        vals = {}
        with io.open(a.estimates) as fh:
            fh.readline()
            for line in fh:
                p = line.split()
                if len(p) >= 2:
                    vals[p[0]] = float(p[1])
        trios = load_trios(a.pedigree)
        raw = trio_reliability(trios, vals)
        print("complete trios used      : %d of %d in the pedigree"
              % (raw["n_trios"], len(trios)))
        print("reliability (raw)        : %.3f +/- %.3f   95%% CI %s"
              % (raw["reliability"], raw["se"], raw["reliability_ci"]))
        print("  total / true / error CV: %.1f%% / %.1f%% / %.1f%%"
              % (raw["total_cv"], raw["true_cv"], raw["error_cv"]))
        print("  excess residual        : %.1f%% CV beyond the Mendelian prediction"
              % raw["excess_resid_cv"])
        if a.pcs:
            adj, r2, npc = adjust_for_pcs(vals, load_pcs(a.pcs, a.n_pc), n_pc=a.n_pc)
            adjr = trio_reliability(trios, adj)
            print("reliability (PC-adjusted): %.3f +/- %.3f   [%d PCs, %.0f%% of "
                  "variance removed]" % (adjr["reliability"], adjr["se"], npc, 100 * r2))
            if adjr["reliability"] < raw["reliability"] - 2 * adjr["se"]:
                print("  NOTE: the raw figure was inflated by error shared within "
                      "families; report the adjusted one.")
        null = permuted_trio_null(trios, vals, n_perm=a.perm)
        print("permuted-family null     : %+.3f +/- %.3f over %d permutations"
              % (null["mean"], null["sd"], null["n_perm"]))
    elif a.cmd == "selftest":
        sys.exit(0 if selftest() else 1)


if __name__ == "__main__":
    main()
