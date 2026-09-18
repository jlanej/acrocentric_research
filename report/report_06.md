# Part IV · Challenges, proposals and an agenda
%SUB What is unresolved, what we propose doing about it, and what we think is true but unproven

## 20. Twelve unresolved challenges

These are ordered by how much we think solving them would unlock, not by difficulty.

**20.1 There is no portable coordinate system.** Every short-arm coordinate is a statement about one assembly. Because array lengths differ between haplotypes by megabases, coordinates do not transfer between individuals even in principle. This blocks cross-study comparison, database curation, clinical reporting and meta-analysis simultaneously. Addressed in 22.1.

**20.2 There is no truth set.** There is no GIAB-style benchmark for the acrocentric short arms — no truth variants, no confident regions, no stratifications. Consequently no caller has ever been evaluated here, and every performance claim in this region is unverified. The medically-relevant-gene benchmarks deliberately stopped short of these regions [@wagner2022cmrg,@olson2023]. Addressed in 22.7.

**20.3 The rDNA arrays remain unassembled, and no current technology will assemble them.** A 45 kb near-identical unit tiled over up to 3.6 Mb cannot be resolved by 20 kb or 100 kb reads. Morph reconstruction [@rautiainen2024ribotin] gives the unit spectrum but not the array order or an exact count. Solving this requires either reads approaching megabase length, single-molecule physical mapping of whole arrays, or an in-situ approach.

**20.4 Assignability is unmeasured.** We have mappability tracks for the reference. We have no metric for "probability that a read from this individual is assigned to the correct arm and haplotype", and therefore no principled way to set a confidence threshold. Addressed in 22.3.

**20.5 Population-scale short-arm data do not exist.** HPRC and HGSVC assemblies are the best available and are incomplete across the short arms [@liao2023,@logsdon2025hgsvc,@porubsky2023gaps]; the deepest short-arm dataset is a single pedigree of 16 informative samples [@lin2026]. Allele frequencies, array-length distributions, hybrid-arm frequencies and morph frequencies are effectively unknown at population scale. Nothing in this region can be interpreted clinically without that denominator.

**20.6 Satellite annotation is heuristic.** Current practice maps consensus sequences and merges hits with length filters [@lin2026]. There is no probabilistic classifier, no calibrated confidence, and no agreed taxonomy boundary between families (the "diverged HSat3" block on 15p at 93% identity to the main array is a case in point). Different groups will annotate the same assembly differently.

**20.7 Array-level variation has no data model.** There is no standard format, no exchange representation, no database schema for "this haplotype's HSat3 array is 4.2 Mb with this unit spectrum". Addressed in 22.4.

**20.8 Statistical genetics for array traits is undeveloped.** Array lengths are quantitative, highly mutable, non-biallelic and measured with substantial error. Standard GWAS and burden machinery does not apply. Mutation rates high enough to produce recurrent identical-by-state alleles break the assumptions of identity-by-descent methods.

**20.9 The functional link is missing.** We can now measure short-arm genotype. We almost never measure the corresponding phenotype — rRNA output, nucleolar number and morphology, ribosome composition, kinetochore function. Without paired functional readouts, short-arm variation cannot be interpreted.

**20.10 Clinical translation has not begun.** No clinical laboratory reports acrocentric short-arm variation from sequencing. Robertsonian carriers are detected by karyotype, an assay from 1959, because no sequencing-based test has been validated. Addressed in 22.6.

**20.11 Reference panels are built from cell lines.** Satellite arrays and rDNA drift in culture, and the community's reference materials — including much of HPRC — are lymphoblastoid lines. The extent to which current "reference" short-arm haplotypes reflect primary genomes is unknown and should be measured.

**20.12 Cost and compute.** Complete short arms currently require three sequencing technologies per sample plus substantial compute. Until the requirement falls to one technology, population-scale short-arm genomics will not happen.

## 21. Connections we think are real but unproven

The following are inferences and hypotheses, not findings. We state each with the prediction that would test it, so that they can be disposed of efficiently.

:::idea 21.1 Synapsis failure as a general mutational mechanism
The observed coincidence on the acrocentric short arms — no allelic crossover, tenfold elevated SNV rate, mismatch-repair-like and oxidative signatures — has been interpreted as causal: without synapsis there is no homologous template, so meiotic double-strand breaks and mismatches are repaired by lower-fidelity routes [@lin2026]. If that mechanism is general, it predicts a genome-wide inverse relationship between local synapsis competence and de novo mutation rate, independent of sequence composition.

The decomposition in section 5.4 is a first partial test and is consistent with the positional model supplying a floor: non-satellite sequence inside the arm mutates 5.6-fold faster than autosomal euchromatin, which composition alone does not explain.

**Predictions.** (i) Regions heterozygous for large inversions should show locally elevated DNM rates in the heterozygous parent's transmissions, relative to the same region in homozygous parents — this is directly testable in existing long-read pedigrees and controls for sequence composition, because the same sequence is compared. (ii) The Yq12 rate (1.99 × 10<super>-7</super> [@lin2026]) should be explained by the same mechanism rather than by satellite composition alone, and should therefore be reproduced in autosomal satellite blocks that *do* synapse — if HSat3 on 1q, 9q or 16q shows a *lower* DNM rate than acrocentric HSat3, composition is not the driver and synapsis is. This is the cleanest single experiment we can propose from this report, and it requires no new sequencing beyond existing pedigree data re-analysed with a complete reference.
:::

:::idea 21.2 The acrocentric short arms are five autosomal non-recombining haplotype systems
If allelic crossover is genuinely absent across the p-arm, each short arm is inherited as a single linkage unit — an autosomal analogue of the male-specific Y or the mitochondrial genome, but present in two copies in every individual of both sexes, five times over. Combined with the tenfold mutation rate, each arm accumulates roughly 10× more substitutions per generation than euchromatin with no recombination to erase the genealogy.

**Consequences, if true.** Short arms should carry deep, well-resolved, population-structured haplogroup phylogenies; they should be exceptionally informative for recent demography and relatedness; ectopic exchange should appear as discrete reticulation events in those phylogenies rather than as a continuous recombination background; and hybrid arms should be datable.

**Predictions.** (i) Linkage disequilibrium measured across a short arm (using Regime I/II variants only) should approach unity within an arm and drop discontinuously at the centromere. (ii) A phylogeny built from Regime I short-arm variants across diverse individuals should be tree-like with a small number of clear reticulations, not net-like. (iii) The five arms should give five partially concordant, partially independent genealogies within the same individual. All three are testable on existing HPRC and HGSVC assemblies restricted to the arm fraction that does assemble — which is a tractable analysis today, and to our knowledge has not been done.
:::

:::idea 21.3 rDNA dosage as a hidden modifier in rare disease
rDNA copy number is heritable, varies two- to three-fold between individuals, is coupled to gene expression and mitochondrial abundance [@gibbons2014,@gibbons2015], and sits upstream of a pathway whose disruption causes tissue-specific developmental disease [@narla2010]. It is also measurable from short-read data that pediatric rare-disease cohorts already hold.

**Prediction.** rDNA copy number should behave as a modifier of penetrance or severity for phenotypes with a nucleolar or ribosomal component — craniofacial, haematopoietic and growth phenotypes in particular — and as a modifier of the phenotypic consequences of aneuploidy. **This is testable at zero sequencing cost** in any existing cohort with short-read WGS and deep phenotyping: estimate rDNA copy number per individual, test for association with severity within diagnosed groups, and test for a copy-number shift in undiagnosed versus diagnosed cases. We consider this the single most actionable idea in this report for a group holding an existing pediatric WGS cohort.
:::

:::idea 21.4 Interlocus gene-conversion tomography between arms
Fourteen of 103 short-arm de novo SNVs matched an allele present elsewhere in the same repeat family, consistent with interlocus gene conversion [@lin2026]. Conversion events are directional and pairwise, which means a pedigree provides a direct measurement of *conversion flux* between specific arms.

**Prediction.** Flux should be strongly asymmetric and should mirror the paralogy structure in Figure 5 — highest between 13p and 21p and between 14p and 22p, lowest involving 15p. Measuring this requires only that conversion-like DNMs be annotated by donor locus rather than discarded as ambiguous, across enough transmissions. It would yield the first direct map of sequence flow between human non-homologous chromosomes, and it provides an independent estimate of the ectopic exchange rate that is far better powered than counting whole-arm recombinants (one event in 13 transmissions).
:::

:::idea 21.5 Nucleolar function as a selective filter on arm haplotypes
If hybrid arms, rDNA morph composition and array lengths affect NOR competence, then arm haplotypes are not neutral, and their population frequencies should deviate from what mutation and drift alone predict.

**Prediction.** rDNA morph frequencies should show less variance between individuals than unit-level mutation and conversion rates predict — i.e. evidence of stabilising selection on morph composition — and individuals with extreme total rDNA copy number should show compensatory shifts in the active/silent ratio measurable by methylation. Both are testable with `ribotin` morph catalogues plus long-read methylation on a modest panel.
:::

:::idea 21.6 Methylation as an orthogonal, cheap haplotype tag
CDR position relative to the α-satellite array start, SST1 methylation level, and satellite-block boundary methylation are all haplotype-specific, stable enough to be transmitted, and obtained free with long-read data [@gershman2022,@lin2026]. Where sequence identity prevents haplotype discrimination, the epigenetic state may not.

**Prediction.** A methylation-derived fingerprint over the short arm should distinguish homologous arms within an individual at higher accuracy than sequence alone in Regime II/III, and should track transmission. If so, "epigenetic haplotype tags" become a practical tool for phasing satellite arrays — a use of methylation data that is currently not exploited at all.
:::

## 22. Seven concrete proposals, and one of them built

These are specifications rather than speculations: each could be built with current data and tooling.

### 22.1 Anchor-relative coordinates for the short arms

**Problem.** Absolute coordinates do not transfer between haplotypes (20.1).

**Proposal.** Address short-arm positions relative to the nearest stable, single-copy-per-arm anchor, with an explicit compartment label and offset:

`13p:DJ+142,318` · `21p:SST1[3]+1,204` · `15p:HOR(S2C15H1L)-88,150`

The anchor set is small and already defined: the subtelomeric ACRO repeat, the DJ, the rDNA array boundaries, the PJ, each named satellite array boundary, and the active HOR array start. Each anchor is identifiable by sequence in any assembly, so the address is computable without a reference. Within an array, the address is (array identifier, unit index, offset in unit) — which is stable against array length change in a way that base coordinates are not.

**What it buys.** Addresses that mean the same thing in two individuals whose arrays differ in length; clinical reports that remain valid when the reference is updated; a natural key for a database of array-level observations. **What it costs.** A translation layer, and a convention for choosing among multiple candidate anchors. We suggest this is a better use of community effort than any further liftover chain.

### 22.2 An acrocentric pangenome module with an explicit ambiguity model

**Problem.** Current human pangenome graphs inherit the assembly gaps and represent paralogy implicitly (20.5, 7.2).

**Proposal.** Build a dedicated acrocentric graph from the subset of assemblies with validated short arms, with three features current graphs lack: (i) nodes annotated with the *set* of arms they are compatible with, not a single position — an explicit representation of the 22–61% of each arm that has a ≥99% twin elsewhere (Figure 5); (ii) array nodes typed as collapsible with a length distribution attached, rather than expanded to one arbitrary length; (iii) rDNA represented as a morph catalogue node rather than a sequence path.

We suggest thinking of the resulting structure as an **onion**: an outer layer of arm-specific sequence, a middle layer shared between one pair of arms, an inner layer shared among all five, and a core (rDNA, HOR arrays) that is not path-resolvable at all. Analysis should proceed layer by layer, reporting results at the outermost layer at which the data can resolve them — rather than forcing everything to the base level and silently mislabelling the inner layers.

### 22.3 An assignability score to replace MAPQ

**Problem.** MAPQ measures uniqueness against one reference haplotype and is misleading here (11.1, 20.4).

**Proposal.** A per-read posterior over arm-of-origin, computed from satellite "dialect" statistics. The pieces exist: satellite dialects — arm-characteristic variant and *k*-mer composition within a satellite family — were described in the original pangenome analysis of these regions [@guarracino2023]; complete assemblies provide training data; the inference is a standard mixture problem.

Concretely: (i) from validated short arms, build for each (arm × satellite family) a profile of discriminative *k*-mers or PSV allele frequencies; (ii) for each read, compute the likelihood under each profile; (iii) report a posterior over the five arms and a null "unassignable" class, fitted by expectation-maximisation over the whole read set so that the individual's own array-size mixture informs the prior. Output a posterior, not a binary placement.

This changes downstream practice: variant calls carry an arm-assignment posterior; the callable region becomes a probabilistic quantity; and a read that is genuinely ambiguous is recorded as ambiguous rather than placed. It is, we think, the most important missing piece of infrastructure for this region.

### 22.4 An array-level variant representation

**Problem.** VCF cannot represent Regime III variation (13.2, 20.7).

**Proposal.** A companion record type — call it an **array-level feature (ALF)** record — keyed by anchor-relative array identity rather than by base coordinate. Minimum fields:

| Field | Content |
| --- | --- |
| `array_id` | anchor-relative identifier, e.g. `15p:HSat3[1]` |
| `sample`, `haplotype` | sample and haplotype identifier |
| `length`, `length_ci` | estimated array length and interval |
| `units`, `unit_spectrum` | unit count and composition (e.g. HOR *n*-mer frequencies) |
| `method`, `evidence` | estimator, input data type, depth |
| `flags` | collapsed, unresolved, boundary-ambiguous, model-derived |
| `methylation` | mean and profile summary where available |
| `assignability` | posterior that the array is correctly assigned to this arm (22.3) |

Two design decisions matter. First, **length is a measurement with an uncertainty**, not a genotype — the format must carry the interval, and downstream statistics must use it. Second, **the flags are mandatory**, because the commonest error in this field is treating a collapsed array as a short one. An ALF record set is also the natural input to array-trait association testing (20.8).

### 22.5 Per-chromosome rDNA copy number by IGS-barcode deconvolution

**Problem.** Total rDNA copy number is measurable; per-chromosome copy number is not (section 16).

**Proposal.** The intergenic spacer carries far more variation than the transcribed unit, and rDNA arrays are homogenised *within* a chromosome more strongly than between chromosomes [@hori2021,@agrawal2018]. Chromosome-characteristic IGS haplotypes should therefore exist. The proposal is a two-stage method: (i) on a long-read reference panel, use `ribotin` morphs together with DJ/PJ-anchored reads to establish which IGS morphs occur on which arm, building an arm × morph matrix; (ii) in a short-read sample, count morph-discriminative *k*-mers and solve the constrained non-negative deconvolution for per-arm copy number, with total copy number as a constraint.

This is a tractable, well-posed inverse problem, and it would turn per-chromosome rDNA dosage — currently a FISH-only measurement — into something obtainable from existing cohort WGS. The main risk is that arm-specific morphs are not sufficiently distinct in all individuals; that risk is itself worth measuring, and stage (i) measures it.

### 22.6 A sequencing-based Robertsonian carrier test, and an acrocentric liability score

**Problem.** ROB carriers are found by karyotype and missed by sequencing (6.1, 20.10).

**Proposal, part one — detection.** A ROB carrier has lost most of two short arms. Three independent short-read signatures follow: reduced depth over arm-specific satellite *k*-mers for the two involved arms; a reduced total rDNA copy number (two arrays lost); and an altered ratio of arm-specific centromeric HOR *k*-mers. A classifier on these features, trained on known carriers — which cytogenetics laboratories hold in quantity — would give a validated, essentially free screen applicable retrospectively to every short-read cohort in existence. Given a population carrier frequency near 1 in 800, any cohort of 10,000 contains of order a dozen carriers, currently unrecognised.

**Proposal, part two — risk stratification.** Since ROB formation is NAHR between high-identity substrates shared between arms [@delima2025rob,@gerton2024], an individual's risk should depend on their *own* haplotypes' shared-substrate content. Define an **acrocentric liability score** for each arm pair as a sum over shared blocks of (block length) × (identity), computed from that individual's assemblies — essentially a personalised version of Figure 5. The prediction is that the score is elevated in ROB carriers' transmitting parents and in families with recurrent ROBs, and that the population distribution of the score predicts the observed 13;14 > 14;21 > 13;21 hierarchy. This is testable on existing ROB cohorts with long-read sequencing of parents, and it would be the first mechanistic risk model for the commonest human structural rearrangement.

### 22.7 An acrocentric benchmark

**Problem.** No truth set exists (20.2).

**Proposal.** Build one on the CEPH1463 pedigree, for which validated short-arm assemblies, three data types and transmission-based validation already exist [@lin2026,@porubsky2025dnm]. Components:

1. **Confident regions**, defined per haplotype by the callable-region intersection in 14.1, published as anchor-relative intervals and as coordinates on each assembly.
2. **Truth variants** in Regimes I and II, validated by orthogonal technology and by transmission.
3. **Truth array-level features** (ALF records, 22.4) for every array that assembles cleanly, with lengths and unit spectra.
4. **Stratifications** that matter here: by regime; by satellite family; by cross-arm paralogy level (using the ≥90/95/99% strata of Figure 5); by distance from the nearest single-copy anchor; and distal versus proximal.
5. **Assignability-aware metrics.** Conventional precision and recall are insufficient: a call that is correct in sequence but assigned to the wrong arm should be scored as a distinct error class, because that is the dominant error mode. We propose reporting *arm-assignment accuracy* alongside genotype concordance, and reporting both stratified by paralogy level.

Even components 1 and 4 alone, published without truth variants, would let the field compare pipelines for the first time.

### 22.8 Implemented: a learned slot grammar, and the placement of unlabelled sequence

Sections 22.1, 22.3 and 22.4 proposed anchor-relative coordinates, an assignability
score and an array-level representation. This section reports a working
implementation of the first and third, an evaluation of the second, and one
instructive failure. Code and tables accompany this report; methods in section 26.

#### The dataset

Per-haplotype satellite annotations were assembled for **944 acrocentric short arms**:
789 from 94 unrelated individuals of HPRC release 2, annotated by the CenSat pipeline,
and 155 from the CEPH1463 pedigree, annotated by the targeted track that additionally
labels ACRO, the DJ, the PJ and SST1. Figure 9 draws all of them.

!FIG fig9_all_arms.png | **Figure 9. Nine hundred and forty-four acrocentric short arms.** One row per haplotype, grouped by chromosome and sorted by length, drawn on each arm's own assembled coordinates. Arms begin at the start of assembled short-arm sequence, which for most HPRC haplotypes is the rDNA array itself because distal sequence is rarely assembled — hence the red left margin. The banding within each chromosome is the conserved block order of section 5.5; the diagonal wedges are size variation in a single array, most visibly the proximal HSat1A block of 13p and 21p and the HSat3 block of 15p. Magenta marks satellite that the annotation calls satellite without naming a family: 120 Mb of it, in 56% of arms. | 0.86

#### A grammar with slots

Because the block order is conserved, an arm can be written as a chain of **slots**,
each carrying a distribution over block families, with alignment allowing a slot to be
skipped and a block to be inserted between slots — the structure of a profile
hidden Markov model, with satellite families in place of residues. Two design choices
make it work here. Arms are split at the rDNA array, the one landmark present on every
arm, and the distal and proximal segments are modelled separately, since the distal
grammar is shared across chromosomes and the proximal grammar is not. And the symbol
for *satellite of unspecified family* emits with probability one from every slot, so
unlabelled blocks are positioned by their context alone and contribute nothing to the
slot distributions they are then compared against.

Slot chains are seeded from an observed arm — the medoid of the group by n-gram
similarity — and refined by Viterbi EM, with slots seeded on a single-copy landmark
held fixed so that the anchors cannot be re-purposed. Fitted models run from 10 slots
(15p proximal) to 21 (14p proximal). Figure 10 shows every arm mapped into the
resulting coordinates.

!FIG fig10_slot_space.png | **Figure 10. The grammar, and every arm expressed in it.** (a) Each row is an arm, each column a slot of that chromosome's fitted proximal grammar, coloured by the family assigned by Viterbi alignment; white means the slot is absent from that arm. Rows are clustered by slot occupancy. Solid columns are slots whose content is effectively fixed across the population; magenta columns are slots consistently occupied by satellite the annotation does not name. (b) The fitted emission distribution of each slot of the 22p proximal grammar — what the grammar permits at each position.

#### Recovering a block's identity from its position

The natural test of a grammar that claims to know what belongs where is to hide a
block and ask it. Masking one block at a time in arms from **held-out individuals**,
the grammar names 89% of masked blocks and is correct for 92.1% of those — 82.3%
of all masked blocks, against 69.7% for a position-bin baseline, 64.1% for majority
class and 35.0% for copying the nearest labelled neighbour (Figure 11a). Trained on
HPRC and tested on the pedigree it holds up: 80.5% overall. The posterior is usable as
a confidence: above 0.95 accuracy is 95–99% (Figure 11b).

!FIG fig11_benchmarks.png | **Figure 11. Evaluation.** (a) Recovery of a masked block's satellite family, for held-out individuals and across datasets, against three baselines. (b) Accuracy by posterior bin; the annotation above each point is the number of blocks in that bin. (c) Assignment of an arm to its chromosome from satellite syntax alone, as a percentage of each true class. (d) Which features carry that signal.

#### Arm-of-origin from syntax alone

Representing each arm as size-weighted counts of ordered block n-grams and classifying
by nearest neighbours, **91.0% of arms from held-out individuals are assigned to the
correct chromosome**, against 21.5% for the largest class (Figure 11c–d). Arm length
alone gives 29.7% and block count 32.3%, so the signal is in composition and order,
not size; order alone, with block sizes discarded, still gives 85.3%.

The errors are the interesting part. They fall almost entirely on 13p–21p (15% of 13p
arms assigned to 21p, 9% the reverse) and 14p–22p (9% and 3%) — the two pairs that
share an active centromeric HOR family, carry the highest cross-arm paralogy
(Figure 5) and account for the recurrent Robertsonian translocations. A classifier
built only on satellite block syntax rediscovers the pairing structure of the
acrocentric community without being told about it.

#### Placing the unlabelled 120 Mb

Applied to the 735 blocks that carry no family assignment, the grammar gives 92% of
them a **reproducible slot address** — a position in the canonical arm, which is the
anchor-relative coordinate of section 22.1 obtained by fitting rather than by decree —
and names 88% of them. The results divide into three honest categories.

**Validated.** 319 blocks totalling 55.1 Mb, median 147 kb, are called ACRO repeat.
This is a vocabulary gap rather than a discovery: the CenSat annotation has no ACRO
label, so genuine subtelomeric ACRO sequence in HPRC assemblies necessarily lands in
the unspecified bin. Masking *labelled* ACRO blocks in pedigree arms recovers them at
85.5% (n=119), which is the evidence that these calls are real.

**Plausible but not validated.** 98 blocks on 21p, median 71 kb, are called SST1 —
squarely inside the published 21–104 kb range for SST1 arrays, at a slot the pedigree
annotation labels SST1, on a chromosome known to carry it. But the corresponding
masking test fails: masked pedigree SST1 blocks are recovered at 0% because pedigree
and HPRC arms take different Viterbi paths through the grammar. The call is
well-motivated and independently consistent in size, and it is not demonstrated. It
would be settled in an afternoon by aligning the SST1 consensus to those coordinates.

**Recurrent and unnamed.** Six slots hold unspecified satellite in at least a quarter
of their chromosome's arms (`data/unknown_slot_summary.csv`). The most prominent is a
proximal slot on 13p present in 62% of arms, 8.4 Mb in total, median block size 85 kb,
for which the grammar offers no confident family (median posterior 0.49); 14p has a
counterpart at 38.5% of arms and median 75 kb. Given that SST1 is documented on 13p,
14p and 21p, and that the equivalent 21p slot *is* confidently called SST1 with a
matching size distribution, the parsimonious reading is that these are the SST1 arrays
of 13p and 14p, unnamed because the annotation vocabulary omits them. The alternative
— that they are a distinct, positionally fixed, unnamed satellite — is the more
interesting possibility and is equally testable.

:::note What this does and does not deliver
It delivers a fitted, inspectable description of where each block type belongs on each
acrocentric short arm; a coordinate for any block, including unlabelled ones, that is
defined by grammatical position rather than by base offset in one assembly; a
calibrated family prediction for unlabelled satellite where the grammar has evidence;
and a chromosome assignment from syntax alone whose errors are biologically
interpretable.

It does not deliver a sequence-level method. Everything here operates on annotation
tracks, so it inherits their errors and their vocabulary, and it cannot discover a
satellite family that no annotation names. The obvious next step is to fit the same
grammar to blocks defined directly from sequence — by *k*-mer composition of each
interval rather than by a label — which would make the "recurrent and unnamed"
category resolvable instead of merely catalogued.
:::

#### The failure worth reporting

**The representation does not transfer between annotation pipelines.** Against 91.0% within a single
pipeline, chromosome assignment drops to 40.8% training on the pedigree track and
testing on HPRC's CenSat track — still above the 21.5% largest-class rate but far below
usable — and to 23.2%, indistinguishable from chance, in the reverse direction. Masked-block
recovery degrades in the same direction (33.8% training on the pedigree, against 80.5%
training on HPRC). The cause is granularity: the two pipelines cut the same sequence
into different numbers of blocks, so n-gram profiles and Viterbi paths are not
comparable across them, and the pedigree's 155 related arms are in any case thin
training data.

This is a concrete, measured obstacle to the assignability score proposed in section
22.3, and it sharpens that proposal. A portable arm-assignment tool cannot be built on
annotation output as currently distributed; it needs either a harmonised satellite
annotation applied uniformly to every assembly, or features computed from sequence
rather than from labels. Of the two, the second is more likely to work, and it is the
direction we would now recommend.

## 23. Pedigree and personal parental references

Every method in sections 9 to 19 assumes the analyst has one genome and a reference built
from strangers. A family breaks that assumption, and it breaks it harder on the
acrocentric short arms than anywhere else in the genome. This section sets out why,
measures how far it goes, and proposes what to build on it. Two of the measurements are
new to this report.

### 23.1 Why a family is worth more here than anywhere else

**The arms do not recombine, so they are inherited as objects.** Across 107 informative
transmissions in the CEPH1463 pedigree, not one allelic crossover was found in the 66 Mb
of short-arm sequence; the 19 events detected in the same analysis all fell in the
proximal q-arms, 18 of 19 of them, with a 2.3-fold enrichment within 5 Mb of the
centromere [@lin2026]. Zero events in 107 transmissions puts the 95% upper bound on the
per-transmission crossover probability at 2.8%, which is an upper bound of 0.042 cM/Mb
against a sex-averaged genome average a little above 1 cM/Mb [@palsson2025] - at least
25-fold lower, and consistent with zero.

The analytical consequence is not "recombination is suppressed", which was already
known. It is that **a parental short arm is transmitted to the child as a unit.**
Everywhere else in the genome a parental haplotype is cut into pieces every generation,
which is why a personal reference decays immediately and why population references are
built at all. In sequence that does not cross over, the parent's arm *is* the child's
arm, altered only by de novo change.

**How much altered is a number we can put down.** At the measured short-arm de novo SNV
rate of 1.33 × 10<super>-7</super> per bp per generation [@lin2026], all five arms
together - 66.07 Mb - are expected to acquire **8.8 substitutions per transmission**,
plus about 0.6 structural events (8 de novo SVs were observed across 13 transmissions,
all in satellite). A transmitting parent's assembled arm therefore differs from the
child's arm at roughly 1.3 × 10<super>-7</super> per bp, where a population reference
differs from any individual at of order 10<super>-3</super> per bp in unique sequence -
a factor of about 7,500, before counting the 83% of the cytogenetic short arms that
GRCh38 does not represent at all (Figure 12d).

:::note
This is the general principle, and it is worth stating in its own right: **personal
references are only durable in sequence that does not recombine.** Elsewhere they are a
one-generation convenience. On the acrocentric short arms, in Yq12, and in the male-specific
region of the Y, a parental assembly is a reference with a decay rate of 10<super>-7</super>
per bp per generation. These are exactly the regions where population references work
worst, which is a fortunate coincidence and an under-exploited one.
:::

### 23.2 Measured: the parts list is transmitted intact

Section 5.5 established that block order is conserved across unrelated individuals while
block sizes vary up to 9.5-fold. That leaves the question a family answers directly: when
a child's arm differs in size from the population, does it match the parent it came from?

We tested this on the pedigree's own annotation - 156 arms with block structure, of which
the eight assembled children of NA12877 and NA12878 contribute 70 arms that can be
compared with both parents. Each child arm was aligned against the candidate arms of each
parent, at most two per parent, so the comparison takes a minimum over the same number of
candidates on each side. Because NA12877 and NA12878 are unrelated, the non-transmitting
parent is an internal population control. Arms were anchored at the proximal end, which is
the end that assembles, using a semi-global alignment with free terminal gaps distally
(methods in section 26). Two quantities are reported per pair: **block-order identity**,
the fraction of aligned block positions with matching family, and **array length ratio**,
the median fold difference in length over matching positions.

!FIG fig12_transmission.png | **Figure 12. What a parent's arm is worth as a reference for a child's.** (a) Block-order identity between a child's arm and the closest candidate arm of the transmitting parent, the non-transmitting parent, and an unrelated individual. Child arms n = 70 from 8 children (paired by arm; Wilcoxon signed-rank); unrelated pairs n = 8,793 from HPRC release 2, same chromosome, order-only. (b) Median array length ratio over matching block positions, same 70 arms. A ratio of 1.00 means every shared array has the same length in parent and child. (c) Fourth-generation arms scored against the great-grandparental couple. Each of the seven individuals has one married-in parent absent from the data, so half of their arms are expected to have no match in the pedigree - the grey mode is that negative control, not a failure. (d) Expected per-base difference between a sample and its reference, for a population reference in unique sequence and for the transmitting parent's arm; the latter is the de novo rate itself.

**Order is transmitted exactly.** Median block-order identity to the transmitting parent
is **1.000**, against 0.833 for the non-transmitting parent (P = 2 × 10<super>-11</super>,
Wilcoxon signed-rank, n = 70); 58.6% of child arms match a parental arm with no block
difference at all, against 10.0% for the non-transmitting parent. The non-transmitting
parent sits exactly where an unrelated individual sits - 0.833 against a median 0.850
across 8,793 same-chromosome pairs of unrelated HPRC arms - which is the control working
as intended: a parent who did not transmit the arm is worth no more than a stranger.

**Length is transmitted exactly too, and that is the load-bearing result.** The median
array length ratio to the transmitting parent is **1.00**, against 1.27 for the
non-transmitting parent (P = 8 × 10<super>-8</super>). Array length is the most variable
property of these arms between individuals and the one that no reference can supply; it
is inherited without measurable change. A parent's arm predicts not just which blocks the
child has and in what order, but how long each one is.

**Where block resolution runs out.** For 32 chromosomes both child haplotypes are
assembled, and 24 of them (75%) resolve cleanly into one paternal and one maternal arm.
The eight failures are not data-poor - they have *more* aligned blocks than the successes
(median 12 against 8) - and in a quarter of the arms involved both parents' candidates
score identically. These are ties, not errors: at the resolution of a satellite block
annotation, the two parents sometimes carry indistinguishable arms. Phasing them needs
sequence, not labels.

**A positive control from the fourth generation.** Each of the seven fourth-generation
individuals has one parent inside the pedigree and one married in and absent from the
data, so half of their arms should trace to the great-grandparental couple and half
should look unrelated. Observed: **52% (24 of 46)** match a great-grandparental arm at
0.9 identity or better (Figure 12c). The measurement tracks descent, not chromosome
identity - if it were responding to chromosome-level similarity alone, all 46 would have
matched.

**And a caution.** What is measured here is the annotation, not the sequence. Identity of
1.000 at block resolution means the parts list and the array lengths agree; it is
consistent with megabases of internal sequence difference, and it says nothing about unit
composition inside an array. It also cannot identify *which* individual transmitted an
arm: a median of 2 of the 8 siblings are indistinguishable per matched arm, because
siblings inherit the same parental arms. An arm identifies a haplotype, not a person.

### 23.3 For assembly

**What a parental assembly supplies that a reference cannot.** Arm assembly currently
fails in three ways (section 9): arrays collapse, paths through the satellite graph are
ambiguous, and distal sequence is orphaned into unplaced contigs. A parental assembly
addresses all three, because section 23.2 says the answer for the child is the parent's
answer: the number of blocks, their order, and their lengths. This converts arm assembly
from an open de novo problem into a guided one, with the guide known to be right to
within 8.8 substitutions.

**Trio binning is not the mechanism, and should not be assumed to be.** Trio-based
phasing assigns reads by parent-specific *k*-mers [@cheng2024hifiasmul,@antipov2025verkko2],
and parent-specific *k*-mers are exactly what satellite arrays are short of - the same
property that gives 13p a 15.6% unique-anchor fraction at 150 bp (Figure 3) limits the
density of markers that distinguish two parental HSat arrays. Trio mode will help most
where it is least needed. **We recommend that trio-phased assemblies report
haplotype-specific marker density per array, not genome-wide**, so that a reader can see
which arrays were actually phased by markers and which were phased by graph topology.

:::idea
**Pedigree-consistent assembly.** Assemble a family jointly rather than each member
separately, and impose inheritance as a constraint on paths: each child's arm path must
be a parental arm path, plus at most the expected de novo budget (of order 9 substitutions
and 0.6 structural events across all five arms). This is a constraint no single-sample
assembler has access to, and it applies precisely where single-sample assemblers are
least confident - the choice between two comparably supported traversals of a satellite
array. In a four-generation pedigree an arm passes three transmissions, so the same
molecule is assembled three or four times from independent data; that is a replicate
structure no single genome provides, and it should be used to resolve paths rather than
merely to check them afterwards.
:::

### 23.4 For alignment, including what it does not fix

This is where the intuition fails, and the failure is quantitative.

**A personal reference does not improve short-read mappability.** The
minimum-unique-anchor ladder of Figure 3 was computed *within a single complete genome*.
It already assumes a perfect, personal, gap-free reference; every anchor it fails to find
is lost to self-similarity inside that genome, not to divergence from a stranger. The
ladder therefore *is* the personal-reference ceiling: at 150 bp, 15.6% of 13p and 45.4%
of 14p carry a unique anchor, and a parental assembly does not move those numbers by a
single base. Nor does it touch cross-arm paralogy, which puts 22.5% to 60.5% of each arm
within 99% identity of sequence on a different acrocentric (Figure 5) - a read ambiguous
between 13p and 21p in CHM13 is ambiguous between them in your mother too.

What a personal reference removes is the other two failure modes: **absent sequence** -
GRCh38 represents 83% of these arms as N, and even CHM13 supplies one person's array
lengths - and **allelic divergence**, which it reduces by about 7,500-fold. Both matter
enormously for long reads and for assembly, and neither is what limits short reads here.

**So the gain is realised at read lengths that can anchor.** Three practical
consequences:

1. Align long reads to a **personal diploid arm panel** - the two parents' ten assembled
   arms - rather than to a reference. Arm-of-origin then has ten candidates drawn from the
   individual's own family, and the residual ambiguity is the paralogy of section 11
   rather than paralogy plus divergence plus missing sequence. Report it as such.
2. The anchor-relative coordinates proposed in section 22.1 become **family-relative
   coordinates**, which is the only system in which a child's array and the parent's array
   it came from have a common origin. Comparing them through a reference coordinate is
   what makes inherited arrays look variant.
3. Where a child's arm has no parental match at block level, that is a finding rather than
   a coordinate problem: it is a candidate de novo rearrangement, or an assembly error in
   one of the two individuals, and the two can be told apart by read-level evidence
   [@mastoras2025].

### 23.5 For analysis

**De novo calling has a hard prior.** Only about 8.8 short-arm substitutions per
transmission exist to be found. Any pipeline reporting hundreds of de novo variants across
these 66 Mb in a trio is reporting errors, and the expected count is a sharper acceptance
test than any filter threshold. The same arithmetic makes the region attractive rather
than hopeless: the de novo rate is tenfold elevated [@lin2026], so the signal-to-noise of
a *correct* pipeline is better here than in euchromatin, once the errors are controlled by
transmission.

**Array-level Mendelian checks, where no variant caller works.** The measurement in
Figure 12b licenses a quality metric that needs no variant calls at all. Array lengths are
inherited unchanged, so for every array in a child's assembly, compare its length - and
its unit spectrum, where `ribotin`-style morph deconvolution is available
[@rautiainen2024ribotin] - with the transmitting parent's. A change is a real array
mutation or an assembly error, and in current data it is usually the latter. We propose
the **array-level Mendelian violation rate** as a reportable QC statistic for acrocentric
assemblies, computed on the ALF records of section 22.4. It is applicable in Regimes II
and III, where precision and recall against a truth set are undefined, and it is the only
validation statistic we know of that works on a collapsed array.

:::note
**The trap that makes trio concordance insufficient.** Satellite assembly errors are
driven by sequence context - array length, unit divergence, coverage dips at the same
repeat - and the parent and child carry the *same* sequence context, because the arm was
inherited intact. The same collapse therefore happens in both assemblies, produces the
same wrong length, and **passes every Mendelian check**. Transmission-based validation is
blind to precisely the error class that dominates these regions. It must be paired with
read-level evidence in both individuals [@mastoras2025], and a benchmark built on
transmission alone (section 22.7) will over-report its own accuracy. We think this is the
single most likely way for a well-intentioned acrocentric truth set to be wrong.
:::

**Clinical uses that need no new method.** For a Robertsonian carrier, parental assemblies
identify which two arms fused and which sequence was lost, rather than inferring it from
the karyotype [@delima2025rob]. For the imprinted disorders that arise from uniparental
disomy of chromosomes 14 and 15 - Temple, Kagami-Ogata, Prader-Willi and Angelman
syndromes [@nicholls2001] - the parent of origin is the diagnosis, and an arm that is
inherited intact is a direct parent-of-origin marker sitting on the same chromosome as
the imprinted locus. Section 23.2 measured how well that works at block resolution: 75%
of chromosomes phase, with the remainder tied rather than wrong.

**One ascertainment caveat, which applies to everything above.** A pedigree is one set of
arm haplotypes from one ancestry. The transmission results here come from a single
European-ancestry family, and the 944-arm HPRC panel of section 22.8 exists to supply the
population context that a family cannot. A pedigree-derived truth set is a truth set for
those haplotypes; it is not a population truth set, and a caller tuned on it will be tuned
on five arms' worth of structure out of a population that varies 9.5-fold in array length.

### 23.6 What to sequence, for a family

| Individual | Recommended data | Why |
| --- | --- | --- |
| Both parents | HiFi 40x + UL-ONT 20x (>100 kb) | The parental arms are the reference everything else uses; they must assemble. |
| Children | HiFi 30x | Order and length come from the parents; the child's data needs only to confirm which arm and catch de novo change. |
| Grandparents | HiFi 30x, if the question is mutation | Each extra generation adds transmissions without adding new haplotypes, which is how a pedigree buys mutation-rate power per genome sequenced [@porubsky2025dnm]. |
| Any individual | ONT or PacBio methylation, from the same reads | NOR activity and rDNA silencing are the functional readout, and it is free with the sequence. |
| Family | Hi-C on one parent only | Scaffolding a parental arm, if it fails to assemble; the children do not need it. |

Two further points of design. **Sequence the parents deeper than the children**, which is
the opposite of the usual trio design where the proband is the priority - here the
parents' assemblies are the shared infrastructure. And **prefer an extra generation to an
extra sibling** when the question is mutation or array dynamics, because siblings
resample the same four parental arms while a third generation supplies new transmissions
of them.

### 23.7 The bridge from families to cohorts

The obvious objection to everything in this section is that it helps only families that
have been sequenced as families. That objection is answerable, and the answer is the one
place where this section and section 22.8 meet.

Pedigree-phased arms are the only source of *verified* long-range arm haplotypes - not
phased by algorithm but by transmission. That makes them the right training set for
imputation, and Figure 12b says what the imputation target should be. Order is conserved
across the population and therefore carries little information; **length is what varies,
and length is what is inherited exactly**. So the quantity to impute for a cohort
individual is the block-length vector, from short-read *k*-mer content against a panel of
pedigree-phased arms, with the slot grammar of section 22.8 supplying the grammar of
which blocks can occupy which position. Validation is built in: hold out a pedigree
member, impute, and compare with the arm they demonstrably inherited.

**Prediction.** Because array lengths are inherited without change and vary several-fold
between individuals, the block-length vector should be substantially heritable and
imputable from *k*-mer content at short-read depth, while the block *order* will be
nearly uninformative. If that holds, cohort-scale acrocentric analysis does not require
cohort-scale long reads - it requires a few hundred pedigree-phased arms and a panel.
That is a far cheaper experiment than the one the field is currently not doing.

!PAGEBREAK

## 24. Cookbook

### 24.1 If you have short reads only

1. Align to T2T-CHM13v2.0 if starting fresh; otherwise keep GRCh38 and add explicit short-arm exclusions.
2. Carry the accessibility mask; report the 2.82 Mb accessible denominator.
3. Report short-arm SNV/indel calls only within the mask; state that the remainder was not assessed.
4. Compute total rDNA copy number by *k*-mer counting against the KY962518 unit, normalised within cohort.
5. Compute relative satellite array content per family (HSat1A, HSat1B, HSat2/3, β-satellite, SST1, α-satellite HOR by family) by *k*-mer depth.
6. Flag outliers on 4 and 5 for long-read follow-up.
7. Never report absence of variation.

### 24.2 If you have HiFi only (~30×)

1. Assemble with hifiasm (UL mode if any ONT is available) or Verkko2.
2. Expect proximal short-arm sequence to assemble, rDNA to collapse, distal sequence to be orphaned or missing.
3. Anchor contigs through ≥1 Mb of q-arm alignment at least 1 Mb from the centromere.
4. Run NucFreq and Flagger; exclude the union.
5. Annotate satellites, SDs, DJ/PJ, rDNA, telomere.
6. Run `ribotin` for rDNA morphs.
7. Report array lengths with collapse flags, per haplotype.

### 24.3 If you have HiFi + UL-ONT + Hi-C

1. Assemble with Verkko2 using all three data types.
2. Identify pq-scatigs; expect ~97% of haplotypes recovered, ~41% with distal sequence scaffolded.
3. Validate with NucFreq, Flagger, ModDotPlot self-dotplots and `CenMAP` for centromeres.
4. Annotate as in 23.2; add CDR identification and methylation.
5. Produce the callable-region set (14.1) and report it.
6. Compare haplotypes all-vs-all with `minimap2 -x asm20 -c --eqx -D -P --dual=no` and 10 kb sliding-window identity.

### 24.4 If you have a trio or pedigree

Read section 23 first: it is the whole reason a family is worth sequencing here. Then follow 24.3 for every member and the de novo procedure in 14.2. Non-negotiable steps: personalised parental reference; two-technology orthogonal support; the 10 kb sliding-window identity filter against distal/proximal misalignment; manual review; gene-conversion annotation of every candidate. Use any further generation as transmission validation, remembering that shared assembly error passes transmission checks (section 23.5).

### 24.5 If you have a Robertsonian or other rearrangement carrier

1. Verkko2 with Hi-C; the graph, not the alignment, is the primary evidence.
2. Expect a fused contig carrying two long arms; characterise the junction by PSV analysis against both parental or reference arms.
3. Confirm with karyotype and, if available, Strand-seq.
4. Report the breakpoint anchor-relative to SST1 and to the DJ/PJ, not only as reference coordinates.

## 25. Reporting standards for acrocentric analyses

We suggest that any paper making a claim about the acrocentric short arms should report, in methods or a supplementary table, the following minimum set. Most current papers report none of it.

1. Reference assembly *and* version, and annotation versions (CenSat, RepeatMasker, SD call set).
2. Whether rDNA models were masked, and by which intervals.
3. The callable denominator, per sample, per haplotype, per chromosome — and how it was computed.
4. The accessibility or mappability criterion applied, with the track cited.
5. Assembly QC outcome: QV, NucFreq/Flagger-flagged fraction, whether distal sequence was recovered and how it was verified (ACRO repeat, DJ).
6. For array lengths: estimator, collapse-detection method, and uncertainty.
7. For rDNA: unit reference used, masking, estimator, and whether the estimate is total or per-chromosome.
8. Arm-assignment procedure, including the anchoring requirement, and the count of contigs that failed it.
9. For variant calls: the orthogonal support requirement, and the gene-conversion check.
10. Explicit statement of what was *not* assessed.

A single sentence — "the acrocentric short arms were not assessed; *n* Mb excluded" — is a complete and honest fulfilment of this standard for a study that excludes them, and is far preferable to silence.
