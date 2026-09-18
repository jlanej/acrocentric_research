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
fails in three ways (section 9.4): arrays collapse, paths through the satellite graph are
ambiguous, and distal sequence is orphaned into unplaced contigs. A parental assembly
addresses all three, because section 23.2 says the answer for the child is the parent's
answer: the number of blocks, their order, and their lengths. This converts arm assembly
from an open de novo problem into a guided one, with the guide known to be right to
within 8.8 substitutions.

**Trio binning is not the mechanism, and should not be assumed to be.** Trio-based
phasing assigns reads by parent-specific *k*-mers [@cheng2024hifiasmul; @antipov2025verkko2],
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
   individual's own family, and the residual ambiguity is the paralogy of section 11.3
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
