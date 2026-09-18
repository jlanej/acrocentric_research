# acrocentric_research

Characterizing the short arms of the human acrocentric chromosomes — 13p, 14p, 15p, 21p and 22p.

This repository holds a consolidated review of these regions and a set of analyses of
their sequence architecture, built from public data only. Nothing here requires
downloading a genome assembly: every analysis runs off published per-haplotype
annotation tracks, which are small.

## What is here

| Path | Contents |
| --- | --- |
| `report/` | `acrocentric_short_arms_report.pdf` — the review, with methods, figures and references, and the markdown sources it is built from |
| `analysis/armsyntax.py` | library: load arms from two annotation dialects, align them as block strings, cluster them, fit and apply a slot grammar |
| `analysis/run_analysis.py` | end-to-end pipeline reproducing every table in `data/` |
| `analysis/fetch_data.sh` | downloads the public inputs (~105 MB) |
| `analysis/render.py`, `analysis/build_report.py` | the report typesetter |
| `data/` | result tables (CSV/JSON) |
| `figures/` | figures at 330 dpi |
| `docs/DATA.md` | provenance of every input file |

## Headline results

**Reference and callability.** The five short arms total 66.07 Mb in T2T-CHM13v2.0.
GRCh38 represents 83% of that territory as N. Only 4.3% of it lies inside the T2T
short-read accessibility mask, against 86.9% of the autosomal genome. Between 22%
and 61% of each arm has a ≥99%-identical copy on a *different* acrocentric
chromosome. Read length buys most of the arm back — unique-anchor availability rises
from 16–45% at 150 bp to 78–94% at 10 kb — and then plateaus; the residual hard core
is the rDNA array and almost nothing else.

**Mutation.** Normalising the published de novo counts by callable bases per
compartment separates the tenfold short-arm elevation into a position-linked floor
(non-satellite sequence inside the arm, 0.76 × 10⁻⁷ per bp per generation, 5.6× the
autosomal rate) and a further composition-linked component (satellite sequence,
2.03 × 10⁻⁷, 15× autosomal). Rate ratio 2.66 (95% CI 1.76–4.02), *P* = 2 × 10⁻⁶.

**Architecture.** Across 944 haplotype-resolved arms — 789 from 94 unrelated HPRC
release 2 individuals and 155 from one pedigree — block *order* is conserved and
block *size* is not. The landmark backbone (ACRO → distal junction → rDNA →
proximal junction → centromere) has an order consistency of 0.97; β-satellite always
occupies the distal-most satellite slot and HSat1A the slot abutting the centromere,
in 100% of informative arms in both datasets; only the relative order of HSat3,
HSat1B and β-satellite in the mid-proximal region is free. Array sizes vary up to
9.5-fold.

**A slot grammar, and placing unlabelled sequence.** Because the order is conserved,
an arm can be modelled as a chain of slots, each with a distribution over block
families, fitted by Viterbi EM. On held-out individuals the fitted grammar recovers a
masked block's satellite family for 82% of blocks (92% accuracy among the 89% it
names), against 70% for a position baseline and 64% for majority class. Satellite
syntax alone assigns an arm to the correct chromosome 91% of the time (chance 22%),
and its errors fall almost entirely on the 13p↔21p and 14p↔22p pairs — the pairs that
share an active centromeric HOR family. Applied to the 735 blocks that HPRC's
annotation calls satellite without specifying a family, the grammar gives 92% of them
a reproducible slot address and names 88% of them; the ACRO calls are validated at
85% by masking labelled pedigree blocks, while the SST1 calls are corroborated by
block size but not by that test. Six recurrent unspecified-satellite slots are
catalogued in `data/unknown_slot_summary.csv`.

**Transmission.** The arms do not cross over, so a parent transmits one as a unit.
Across 70 arms from eight children of a single couple, block order matches the
transmitting parent exactly in 58.6% of arms (median identity 1.000 vs 0.833 for the
non-transmitting parent, P = 2e-11) and array lengths match with a median ratio of
1.00 vs 1.27 (P = 8e-8) - array length being the most variable property of these arms
between individuals. The non-transmitting parent scores like a stranger (0.833 vs a
median 0.850 over 8,793 unrelated same-chromosome pairs), and 52% of fourth-generation
arms trace to the great-grandparental couple against the 50% expected from one
married-in parent each. A parental arm therefore differs from the child's at the de
novo rate itself, 1.3e-7/bp - about 7,500-fold closer than a population reference. It
does not improve short-read mappability, which is set by self-similarity within any
one genome and is already what the anchor ladder measures.

**What did not work, and matters.** The syntax representation does not transfer
between annotation pipelines. Against 91.0% within a single pipeline (held-out
individuals, chance 21.5%), chromosome assignment drops to 40.8% training on the
pedigree track and testing on HPRC's CenSat track, and to 23.2% - chance - in the
reverse direction. The two pipelines cut the same sequence into different numbers of
blocks, so harmonised satellite annotation, or features computed from sequence rather
than from labels, is a prerequisite for any portable arm-assignment tool.

## Repeat copy number from WGS

`analysis/repeatcn.py` estimates the dosage of multi-copy sequence classes - rDNA 45S
and 5S, the HSat families, alpha/beta/gamma satellite, SST1, ACRO, telomere, and
macrosatellites given a panel - from short-read WGS. Dosage is the one measurement
short reads support well inside these regions, and rDNA copy number has published
associations with blood cell counts and renal function in UK Biobank.

Depth extraction is delegated to mosdepth; the module owns the parts that are easy to
get wrong - the estimator, the denominator, GC correction, technical covariates, and
the controls. Design notes and the class table are in `docs/REPEAT_CN.md`; the target
table is `data/repeat_cn_targets.tsv`.

Batch is handled by residualising on coverage PCs from
[NGS-PCA](https://github.com/jlanej/NGS-PCA) rather than on recorded batch labels. The
PC basis is built from single-copy bins with satellite, segmental duplications and low
mappability excluded, so it is disjoint from every class being measured: the components
can absorb library and bias structure but cannot absorb the dosage of the target.

`python analysis/repeatcn.py selftest` validates the whole thing on simulated counts,
with no data required. It checks copy-number recovery, the shared-denominator artefact
and its removal by a split denominator, GC correction at 58% GC (-72% error to -2%),
sub-region concordance, and PC adjustment of a batch-confounded cohort where the raw
association has the wrong sign (true r = +0.09, raw -0.14, adjusted +0.12).

## Reproducing

```sh
cd analysis
./fetch_data.sh ../work            # ~105 MB of public annotation tracks
python run_analysis.py --work ../work --out ../data --fig ../figures
```

Requires Python 3 with numpy, pandas, scipy and matplotlib. Report rebuild
additionally needs reportlab.

## Caveats

The pedigree arms are related, so conserved order there partly reflects shared
descent; the HPRC replication in unrelated individuals is what carries the fixed
relationships. Both datasets only see arms that assemble, so order variation severe
enough to prevent assembly is invisible. All de novo mutation counts come from a
single pedigree. The CHM13 rDNA arrays are models, and any statistic computed inside
them describes the model rather than human rDNA.

