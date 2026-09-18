"""
armsyntax — represent, align, cluster and model acrocentric short arms as
sequences of satellite blocks ("satellite syntax").

The unit of analysis is an *arm*: an ordered list of (family, start, end)
blocks running from the telomere to the start of the active alpha-satellite
array. Two annotation dialects are supported:

  * "pedigree"  — the Platinum Pedigree AcroMutRecomb `all_samples_censat.bed`
                  track, which carries landmark labels (ACRO, DJ, PJ, SST1)
                  that the standard CenSat vocabulary lacks.
  * "censat"    — the per-haplotype CenSat annotations distributed for HPRC
                  release 2 assemblies.

Both are mapped onto one common alphabet (ALPHABET) so that arms from
different sources are directly comparable.
"""
from __future__ import annotations
import re, glob, os, collections
import numpy as np

# ----------------------------------------------------------------------------
# common alphabet
# ----------------------------------------------------------------------------
#  known satellite families and landmarks, the non-satellite spacer class, and
#  a wildcard for satellite whose family the annotation does not specify.
ALPHABET = ["ACRO", "DJ", "rDNA", "PJ", "SST1", "bSat", "HSat1A", "HSat1B",
            "HSat2", "HSat3", "gSat", "aSat", "NONSAT", "UNK"]
KNOWN_SAT = ["bSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "gSat", "aSat"]
LANDMARK = ["ACRO", "DJ", "rDNA", "PJ", "SST1"]
WILDCARD = "UNK"          # satellite of unspecified family - to be placed
SPACER = "NONSAT"         # non-satellite / transition sequence

# pedigree dialect -> common alphabet
_PED = {
    "ACRO": "ACRO", "DJarm": "DJ", "DJflank": "DJ", "rDNA": "rDNA", "PJ": "PJ",
    "SST1": "SST1", "bSat": "bSat", "HSat1A": "HSat1A", "HSat1B": "HSat1B",
    "HSat2": "HSat2", "HSat3": "HSat3", "aSat": "aSat", "SEQ": "NONSAT",
    # unnamed / generic satellite classes in that track
    "CER": "UNK", "PHR_LOC": "UNK", "yellow1": "UNK", "red2": "UNK",
    "gray1": "UNK", "green2": "UNK", "blue2": "UNK",
    # Y-amplicon targets that cross-match on autosomes: not interpretable here
    "P5-AZFb": None, "P3-AZFc": None, "TEL": None, "GAP": None,
}
# CenSat dialect -> common alphabet
_CENSAT = {
    "bSat": "bSat", "HSat1A": "HSat1A", "HSat1B": "HSat1B", "HSat2": "HSat2",
    "HSat3": "HSat3", "gSat": "gSat", "rDNA": "rDNA",
    "mon": "aSat", "hor": "aSat", "dhor": "aSat", "active_hor": "aSat",
    "mixedAlpha": "UNK", "cenSat": "UNK", "ct": "NONSAT",
    "GAP": None, "TEL": None,
}

ACTIVE_FAMILY = {"S2C13/21H1L": ("chr13", "chr21"),
                 "S2C14/22H1L": ("chr14", "chr22"),
                 "S2C15H1L": ("chr15",)}
CHROM_LEN = {"chr13": 113.6e6, "chr14": 101.2e6, "chr15": 99.8e6,
             "chr21": 45.1e6, "chr22": 51.3e6}


# ----------------------------------------------------------------------------
# block building
# ----------------------------------------------------------------------------
def blocks_from_intervals(ivs, min_bp=50_000, join_bp=50_000):
    """Collapse annotated intervals into non-overlapping blocks.

    ivs : iterable of (family, start, end) on the common alphabet, any order.
    Consecutive intervals of the same family within `join_bp` are merged;
    overlaps between different families are resolved in favour of the block
    already laid down (annotations are near-tiling, so this is rare).
    Blocks shorter than `min_bp` are dropped.
    """
    out = []
    for fam, s, e in sorted(ivs, key=lambda x: (x[1], x[2])):
        if e <= s:
            continue
        if out and fam == out[-1][0] and s <= out[-1][2] + join_bp:
            out[-1][2] = max(out[-1][2], e)
        elif out and s < out[-1][2]:
            if e > out[-1][2]:
                out.append([fam, out[-1][2], e])
        else:
            out.append([fam, s, e])
    return [tuple(b) for b in out if b[2] - b[1] >= min_bp]


class Arm:
    """One acrocentric short arm: ordered blocks from telomere to centromere."""

    __slots__ = ("sample", "hap", "chrom", "source", "blocks", "length",
                 "complete", "tig")

    def __init__(self, sample, hap, chrom, source, blocks, length, tig=""):
        self.sample, self.hap, self.chrom = sample, hap, chrom
        self.source, self.blocks, self.length, self.tig = source, blocks, length, tig
        fams = {b[0] for b in blocks}
        self.complete = ("rDNA" in fams) and (
            "DJ" in fams or "ACRO" in fams or blocks[0][1] < 200_000 if blocks else False)

    @property
    def symbols(self):
        return [b[0] for b in self.blocks]

    @property
    def id(self):
        return f"{self.sample}_h{self.hap}_{self.chrom}"

    def collapsed(self):
        """Block symbols with immediate repeats collapsed."""
        s = self.symbols
        return [x for i, x in enumerate(s) if i == 0 or x != s[i - 1]]

    def __repr__(self):
        return f"<Arm {self.id} {self.length/1e6:.2f}Mb {len(self.blocks)} blocks>"


# ----------------------------------------------------------------------------
# loaders
# ----------------------------------------------------------------------------
def load_pedigree(bed_path, min_bp=50_000):
    """Load arms from the Platinum Pedigree all_samples_censat.bed track."""
    rows = collections.defaultdict(list)
    meta = {}
    for line in open(bed_path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or line.startswith("track"):
            continue
        tig, s, e, lab = f[0], int(f[1]), int(f[2]), f[3]
        m = re.match(r"^([^_]+)_(\d)_haplotype(\d)-(\d+)_(chr\d+)$", tig)
        if not m:
            continue
        meta[tig] = (m.group(1), int(m.group(3)), m.group(5))
        rows[tig].append((lab, s, e))
    arms = []
    for tig, ivs in rows.items():
        sample, hap, chrom = meta[tig]
        asat = [(s, e) for lab, s, e in ivs if lab == "aSat"]
        if not asat:
            continue
        cen_start = max(asat, key=lambda x: x[1] - x[0])[0]
        mapped = []
        for lab, s, e in ivs:
            fam = _PED.get(lab, "UNK")
            if fam is None or s >= cen_start:
                continue
            mapped.append((fam, s, min(e, cen_start)))
        blocks = blocks_from_intervals(mapped, min_bp=min_bp)
        if blocks:
            arms.append(Arm(sample, hap, chrom, "pedigree", blocks, cen_start, tig))
    return arms


def _read_fai(path):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(path)}


def load_censat_haplotype(bed_path, fai_path, min_bp=50_000, len_tol=0.15):
    """Load acrocentric arms from one HPRC-style per-haplotype CenSat BED.

    Scaffolds are assigned to a chromosome by the suprachromosomal-family-2
    active HOR family named in the annotation, disambiguated within a family by
    scaffold length; only chromosome-scale scaffolds are kept, one per
    chromosome, and scaffolds whose centromere lies distally are flipped.
    """
    base = os.path.basename(bed_path)
    sample = base.split("_")[0]
    hap = int(re.search(r"hap(\d)", base).group(1))
    L = _read_fai(fai_path)
    ivs = collections.defaultdict(list)
    active = collections.defaultdict(list)
    for line in open(bed_path):
        if line.startswith("track"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        tig, s, e, lab = f[0], int(f[1]), int(f[2]), f[3]
        ivs[tig].append((lab, s, e))
        if "active_hor" in lab:
            active[tig].append((lab, s, e))
    cands = {}
    for tig, acts in active.items():
        fam = None
        for lab, _, _ in acts:
            m = re.search(r"S2C(13/21|14/22|15)H1L", lab)
            if m:
                fam = "S2C" + m.group(1) + "H1L"
                break
        if fam is None or tig not in L:
            continue
        n = L[tig]
        ch = min(ACTIVE_FAMILY[fam], key=lambda c: abs(n - CHROM_LEN[c]))
        fit = abs(n - CHROM_LEN[ch]) / CHROM_LEN[ch]
        if fit > len_tol:
            continue
        cst, cend = min(a[1] for a in acts), max(a[2] for a in acts)
        if ch not in cands or fit < cands[ch][0]:
            cands[ch] = (fit, tig, n, cst, cend)
    arms = []
    for ch, (_, tig, n, cst, cend) in cands.items():
        rev = (cst / n) > 0.5
        cen_start = (n - cend) if rev else cst
        if cen_start < 1e6:
            continue
        mapped = []
        for lab, s, e in ivs[tig]:
            fam = _CENSAT.get(lab.split("(")[0], "UNK")
            if fam is None:
                continue
            s2, e2 = (n - e, n - s) if rev else (s, e)
            if s2 >= cen_start:
                continue
            mapped.append((fam, s2, min(e2, cen_start)))
        blocks = blocks_from_intervals(mapped, min_bp=min_bp)
        if blocks:
            arms.append(Arm(sample, hap, ch, "censat", blocks, cen_start, tig))
    return arms


def load_censat_dir(bed_dir, fai_dir, min_bp=50_000, limit=None):
    arms = []
    beds = sorted(glob.glob(os.path.join(bed_dir, "*.cenSat.bed")))
    if limit:
        beds = beds[:limit]
    for b in beds:
        stem = os.path.basename(b).split("_hprc")[0]          # SAMPLE_hapN
        fai = glob.glob(os.path.join(fai_dir, stem + "_hprc*.fai"))
        if not fai:
            continue
        try:
            arms += load_censat_haplotype(b, fai[0], min_bp=min_bp)
        except Exception:
            continue
    return arms


# ============================================================================
# 1. syntax alignment: pairwise global alignment of block-symbol strings
# ============================================================================
IDX = {s: i for i, s in enumerate(ALPHABET)}


def default_submat(match_sat=2.0, match_landmark=5.0, match_spacer=0.6,
                   mismatch=-2.5, wildcard=0.0):
    """Score matrix over ALPHABET. Landmarks score high because they are
    single-copy per arm and therefore strong anchors; the wildcard scores
    neutrally against everything so that unspecified satellite neither
    supports nor penalises an alignment."""
    K = len(ALPHABET)
    S = np.full((K, K), mismatch, float)
    for a in ALPHABET:
        for b in ALPHABET:
            i, j = IDX[a], IDX[b]
            if a == WILDCARD or b == WILDCARD:
                S[i, j] = wildcard
            elif a == b:
                S[i, j] = (match_landmark if a in LANDMARK
                           else match_spacer if a == SPACER else match_sat)
    return S


def align(x, y, S, gap_open=-3.0, gap_ext=-0.8, traceback=True):
    """Affine-gap global alignment of two symbol lists. Returns (score, pairs)
    where pairs are (i, j) aligned index pairs. Pure-Python inner loop: the
    matrices are small and this is called O(n^2) times."""
    n, m = len(x), len(y)
    NEG = -1e18
    Srow = [[float(v) for v in row] for row in S]
    xi = [IDX[a] for a in x]; yj = [IDX[b] for b in y]
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG] * (m + 1) for _ in range(n + 1)]
    pM = [[0] * (m + 1) for _ in range(n + 1)]
    pIx = [[0] * (m + 1) for _ in range(n + 1)]
    pIy = [[0] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0.0
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + gap_ext * (i - 1); pIx[i][0] = 1
    for j in range(1, m + 1):
        Iy[0][j] = gap_open + gap_ext * (j - 1); pIy[0][j] = 2
    for i in range(1, n + 1):
        Mi, Mp = M[i], M[i - 1]
        Ixi, Ixp = Ix[i], Ix[i - 1]
        Iyi, Iyp = Iy[i], Iy[i - 1]
        pMi, pIxi, pIyi = pM[i], pIx[i], pIy[i]
        srow = Srow[xi[i - 1]]
        for j in range(1, m + 1):
            s = srow[yj[j - 1]]
            a, b, c = Mp[j - 1], Ixp[j - 1], Iyp[j - 1]
            if a >= b and a >= c:
                Mi[j] = a + s; pMi[j] = 0
            elif b >= c:
                Mi[j] = b + s; pMi[j] = 1
            else:
                Mi[j] = c + s; pMi[j] = 2
            a2, b2 = Mp[j] + gap_open, Ixp[j] + gap_ext
            if a2 >= b2:
                Ixi[j] = a2; pIxi[j] = 0
            else:
                Ixi[j] = b2; pIxi[j] = 1
            a3, b3 = Mi[j - 1] + gap_open, Iyi[j - 1] + gap_ext
            if a3 >= b3:
                Iyi[j] = a3; pIyi[j] = 0
            else:
                Iyi[j] = b3; pIyi[j] = 1
    ends = (M[n][m], Ix[n][m], Iy[n][m])
    state = 0 if ends[0] >= ends[1] and ends[0] >= ends[2] else (1 if ends[1] >= ends[2] else 2)
    score = ends[state]
    if not traceback:
        return score, None
    i, j, pairs = n, m, []
    while i > 0 or j > 0:
        if state == 0:
            pairs.append((i - 1, j - 1)); k = pM[i][j]; i -= 1; j -= 1; state = k
        elif state == 1:
            k = pIx[i][j]; i -= 1; state = 0 if k == 0 else 1
        else:
            k = pIy[i][j]; j -= 1; state = 0 if k == 0 else 2
        if i == 0 and j == 0:
            break
    return score, pairs[::-1]


# ---------------------------------------------------------------------------
# fast order-aware representation: size-weighted block n-grams
# ---------------------------------------------------------------------------
def ngram_features(arms, orders=(1, 2, 3), weight="sqrt_bp", collapse=True,
                   drop_wildcard=False):
    """Vectorise arms as size-weighted counts of ordered block n-grams.

    Each arm contributes, for every window of n consecutive blocks, one count
    against the tuple of their families. Counts are weighted by block size
    (`sqrt_bp`, `bp` or `none`) so that a megabase array is not equivalent to
    a 60 kb one, then L2-normalised. Cosine distance in this space is an
    order-aware, truncation-tolerant arm distance computable for thousands of
    arms at once.
    """
    feats = []
    vocab = {}
    for a in arms:
        blocks = a.blocks
        if collapse:
            merged = []
            for fam, s, e in blocks:
                if merged and merged[-1][0] == fam:
                    merged[-1][2] = e
                else:
                    merged.append([fam, s, e])
            blocks = [tuple(b) for b in merged]
        if drop_wildcard:
            blocks = [b for b in blocks if b[0] != WILDCARD]
        d = collections.Counter()
        for n in orders:
            for i in range(len(blocks) - n + 1):
                win = blocks[i:i + n]
                key = (n,) + tuple(b[0] for b in win)
                bp = sum(b[2] - b[1] for b in win) / n
                w = {"bp": bp, "sqrt_bp": np.sqrt(bp), "none": 1.0}[weight]
                d[key] += w
                if key not in vocab:
                    vocab[key] = len(vocab)
        feats.append(d)
    X = np.zeros((len(arms), len(vocab)))
    for r, d in enumerate(feats):
        for k, v in d.items():
            X[r, vocab[k]] = v
    nrm = np.linalg.norm(X, axis=1, keepdims=True)
    X = X / np.where(nrm == 0, 1, nrm)
    return X, vocab


def cosine_distance(X):
    G = X @ X.T
    np.clip(G, -1, 1, out=G)
    return 1.0 - G


def syntax_distance(arms, S=None, collapse=True, **kw):
    """Symmetric normalised alignment distance between every pair of arms.

    d(a,b) = 1 - 2*s(a,b) / (s(a,a) + s(b,b)),  clipped to [0, 1].
    """
    S = default_submat() if S is None else S
    seqs = [a.collapsed() if collapse else a.symbols for a in arms]
    n = len(seqs)
    self_s = np.array([align(s, s, S, **kw)[0] for s in seqs])
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sij = align(seqs[i], seqs[j], S, **kw)[0]
            denom = self_s[i] + self_s[j]
            d = 1.0 - (2.0 * sij / denom) if denom > 0 else 1.0
            D[i, j] = D[j, i] = min(max(d, 0.0), 1.0)
    return D


def reestimate_submat(arms, S=None, collapse=True, pseudo=1.0, **kw):
    """Data-driven substitution matrix: log-odds of two families occupying
    aligned positions, in the style of a BLOSUM, estimated from one round of
    all-pairs alignment under a seed matrix."""
    S0 = default_submat() if S is None else S
    seqs = [a.collapsed() if collapse else a.symbols for a in arms]
    K = len(ALPHABET)
    obs = np.full((K, K), pseudo)
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            for (a, b) in align(seqs[i], seqs[j], S0, **kw)[1]:
                u, v = IDX[seqs[i][a]], IDX[seqs[j][b]]
                obs[u, v] += 1; obs[v, u] += 1
    p = obs / obs.sum()
    marg = p.sum(1)
    exp = np.outer(marg, marg)
    return np.log2(p / exp)


# ============================================================================
# 2. slot grammar: a profile of ordered slots, fitted by Viterbi EM
# ============================================================================
SEED_SLOTS = ["NONSAT", "bSat", "HSat3", "HSat1A", "HSat3", "NONSAT", "ACRO",
              "NONSAT", "DJ", "rDNA", "PJ", "NONSAT", "bSat", "HSat3",
              "HSat1B", "HSat3", "bSat", "HSat3", "SST1", "HSat1A", "aSat"]
EMIT_SYMBOLS = [s for s in ALPHABET if s != WILDCARD]


def collapse_blocks(blocks):
    out = []
    for fam, s, e in blocks:
        if out and out[-1][0] == fam:
            out[-1][2] = e
        else:
            out.append([fam, s, e])
    return [tuple(x) for x in out]


def split_segments(arm, drop_terminal_asat=True):
    """Split an arm at its rDNA array into (distal, proximal) block lists.

    Returns None when the arm has no rDNA anchor. The rDNA array is the one
    landmark present on every acrocentric short arm, so splitting on it gives
    two segments that can be modelled independently - which is what the
    observed grammar requires, the distal order being shared across
    chromosomes and the proximal order being chromosome-specific.
    """
    rd = [b for b in arm.blocks if b[0] == "rDNA"]
    if not rd:
        return None
    r0 = min(b[1] for b in rd); r1 = max(b[2] for b in rd)
    distal = collapse_blocks([b for b in arm.blocks if b[2] <= r0])
    prox = [b for b in arm.blocks if b[1] >= r1]
    if drop_terminal_asat:
        prox = [b for b in prox if b[0] != "aSat"]
    return distal, collapse_blocks(prox)


def seed_from_medoid(seqs, lo=0.60, hi=0.92):
    """Choose a seed slot chain: among sequences whose length falls between the
    `lo` and `hi` quantiles, take the one with the smallest mean n-gram cosine
    distance to all others. This is a seed-and-refine strategy - the seed is an
    observed arm rather than a hand-written consensus."""
    seqs = [s for s in seqs if s]
    if not seqs:
        return []
    lens = np.array([len(s) for s in seqs])
    a, b = np.quantile(lens, lo), np.quantile(lens, hi)
    cand = [i for i, s in enumerate(seqs) if a <= len(s) <= b] or list(range(len(seqs)))

    def prof(s):
        d = collections.Counter()
        for n in (1, 2):
            for i in range(len(s) - n + 1):
                d[(n,) + tuple(s[i:i + n])] += 1.0
        return d
    P = [prof(s) for s in seqs]
    keys = {k for p in P for k in p}
    ki = {k: i for i, k in enumerate(keys)}
    X = np.zeros((len(seqs), len(ki)))
    for r, p in enumerate(P):
        for k, v in p.items():
            X[r, ki[k]] = v
    X /= np.where(np.linalg.norm(X, axis=1, keepdims=True) == 0, 1,
                  np.linalg.norm(X, axis=1, keepdims=True))
    D = 1 - X @ X.T
    best = min(cand, key=lambda i: D[i].mean())
    return list(seqs[best])


class SlotGrammar:
    """A left-to-right profile over the block alphabet.

    Slots are positions in the canonical arm; each slot carries a
    distribution over block families. An arm is aligned to the profile by
    Viterbi, allowing a slot to be skipped (deletion) and a block to be
    inserted between slots. The wildcard symbol emits with probability 1 from
    every slot, so unspecified satellite is placed by its context alone and
    the slot's distribution over *known* families becomes a prediction of what
    it is.
    """

    def __init__(self, seed=None, conc=6.0, pseudo=0.25,
                 t_del=-2.2, t_ins=-2.6, t_ins_ext=-0.7):
        self.slots = list(seed or SEED_SLOTS)
        self.L = len(self.slots)
        self.syms = list(EMIT_SYMBOLS)
        self.si = {s: i for i, s in enumerate(self.syms)}
        self.t_del, self.t_ins, self.t_ins_ext = t_del, t_ins, t_ins_ext
        E = np.full((self.L, len(self.syms)), pseudo)
        for j, f in enumerate(self.slots):
            if f in self.si:
                E[j, self.si[f]] += conc
        self.E = E / E.sum(1, keepdims=True)
        self.bg = np.full(len(self.syms), 1.0 / len(self.syms))

    # -- emission log-probability, with the wildcard as a free pass ----------
    def _logE(self, j, sym):
        if sym == WILDCARD:
            return 0.0
        return float(np.log(self.E[j, self.si[sym]] + 1e-12))

    def _logbg(self, sym):
        if sym == WILDCARD:
            return 0.0
        return float(np.log(self.bg[self.si[sym]] + 1e-12))

    def viterbi(self, syms):
        """Align one arm to the profile. Returns (score, assignment) where
        assignment[i] is the slot index for block i, or -1 for an insertion."""
        n, L = len(syms), self.L
        NEG = -1e9
        # V[i][j]: best score having consumed i blocks and j slots, block i-1
        # matched to slot j-1 (state M) or inserted (state I)
        VM = np.full((n + 1, L + 1), NEG)
        VI = np.full((n + 1, L + 1), NEG)
        BM = np.zeros((n + 1, L + 1), np.int8)
        BI = np.zeros((n + 1, L + 1), np.int8)
        VM[0, 0] = 0.0
        for j in range(1, L + 1):                       # leading deletions
            VM[0, j] = VM[0, j - 1] + self.t_del
            BM[0, j] = 3
        for i in range(1, n + 1):
            for j in range(0, L + 1):
                # insertion: consume block i, stay at slot j
                opts = (VM[i - 1, j] + self.t_ins, VI[i - 1, j] + self.t_ins_ext)
                k = int(np.argmax(opts))
                VI[i, j] = opts[k] + self._logbg(syms[i - 1]); BI[i, j] = k
                if j == 0:
                    continue
                e = self._logE(j - 1, syms[i - 1])
                cand = (VM[i - 1, j - 1] + e,            # 0 match from match
                        VI[i - 1, j - 1] + e,            # 1 match from insert
                        VM[i, j - 1] + self.t_del,       # 3 slot deleted
                        VI[i, j - 1] + self.t_del)
                kk = int(np.argmax(cand))
                VM[i, j] = cand[kk]
                BM[i, j] = (0, 1, 3, 4)[kk]
        # terminate in whichever state is best at (n, L)
        state = 0 if VM[n, L] >= VI[n, L] else 1
        score = max(VM[n, L], VI[n, L])
        assign = [-1] * n
        i, j = n, L
        while i > 0 or j > 0:
            if state == 0:
                b = BM[i, j]
                if b in (0, 1):
                    assign[i - 1] = j - 1; i -= 1; j -= 1; state = 0 if b == 0 else 1
                elif b in (3, 4):
                    j -= 1; state = 0 if b == 3 else 1
                else:
                    break
            else:
                b = BI[i, j]
                assign[i - 1] = -1; i -= 1; state = 0 if b == 0 else 1
            if i <= 0 and j <= 0:
                break
        return score, assign

    def fit(self, seqs, rounds=8, pseudo=0.25, verbose=False, pin_landmarks=True):
        """Viterbi (hard) EM: align all sequences, re-estimate slot emissions
        from the blocks assigned to each slot, repeat. The wildcard contributes
        no emission counts, so slot distributions describe known families only.

        `seqs` is a list of symbol lists (use `Arm.symbols`, or one segment of
        `split_segments`). With `pin_landmarks`, slots seeded on a single-copy
        landmark keep their seeded emission distribution, which stops EM from
        re-purposing the anchors that make the alignment interpretable.
        """
        frozen = [j for j, f in enumerate(self.slots)
                  if pin_landmarks and f in LANDMARK]
        E0 = self.E.copy()
        hist = []
        for r in range(rounds):
            counts = np.full_like(self.E, pseudo)
            ins = np.full(len(self.syms), pseudo)
            tot = 0.0
            for syms in seqs:
                if not syms:
                    continue
                sc, asg = self.viterbi(syms)
                tot += sc
                for sym, j in zip(syms, asg):
                    if sym == WILDCARD:
                        continue
                    if j >= 0:
                        counts[j, self.si[sym]] += 1
                    else:
                        ins[self.si[sym]] += 1
            self.E = counts / counts.sum(1, keepdims=True)
            for j in frozen:
                self.E[j] = E0[j]
            self.bg = ins / ins.sum()
            hist.append(tot / max(len(seqs), 1))
            if verbose:
                print(f"  round {r+1}: mean score {hist[-1]:.3f}")
            if len(hist) > 1 and abs(hist[-1] - hist[-2]) < 1e-4:
                break
        self.history = hist
        return self

    # -- use -----------------------------------------------------------------
    def slot_table(self, top=3):
        """Human-readable description of each fitted slot."""
        rows = []
        for j in range(self.L):
            o = np.argsort(-self.E[j])[:top]
            rows.append(dict(slot=j, seed=self.slots[j],
                             top=", ".join(f"{self.syms[k]} {self.E[j,k]:.2f}" for k in o)))
        return rows

    def place_unknown(self, arm, restrict=None):
        """Assign every block of an arm to a slot and, for wildcard blocks,
        return the posterior over known satellite families implied by the slot.
        """
        restrict = restrict or KNOWN_SAT
        _, asg = self.viterbi(arm.symbols)
        out = []
        for (sym, s, e), j in zip(arm.blocks, asg):
            rec = dict(arm=arm.id, chrom=arm.chrom, source=arm.source,
                       sym=sym, start=s, end=e, bp=e - s, slot=j)
            if j >= 0:
                p = np.array([self.E[j, self.si[f]] for f in restrict])
                p = p / p.sum() if p.sum() > 0 else p
                k = int(np.argmax(p))
                rec.update(pred=restrict[k], post=float(p[k]),
                           posterior={f: float(v) for f, v in zip(restrict, p)})
            else:
                rec.update(pred=None, post=np.nan, posterior=None)
            out.append(rec)
        return out

    def score_arm(self, arm, n_null=20, rng=None):
        """Log-odds of an arm against the profile, relative to permutations of
        its own block order. Low values flag unusual architecture."""
        rng = rng or np.random.default_rng(0)
        s, _ = self.viterbi(arm.symbols)
        null = []
        sy = list(arm.symbols)
        for _ in range(n_null):
            rng.shuffle(sy)
            null.append(self.viterbi(sy)[0])
        null = np.array(null)
        return s, float(np.mean(null)), float((s - null.mean()) / (null.std() + 1e-9))

