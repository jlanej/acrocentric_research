### 22.8 Implemented: a learned slot grammar, and the placement of unlabelled sequence

Sections 22.1, 22.3 and 22.4 proposed anchor-relative coordinates, an assignability
score and an array-level representation. This section reports a working
implementation of the first and third, an evaluation of the second, and one
instructive failure. Code and tables accompany this report; methods in section 25.

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

The errors are the interesting part. They fall almost entirely on 13p↔21p (15% of 13p
arms assigned to 21p, 9% the reverse) and 14p↔22p (9% and 3%) — the two pairs that
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

**The representation does not transfer between annotation pipelines.** Trained on the
pedigree track and tested on HPRC's CenSat track, chromosome assignment collapses from
91.0% to 40.8%, and in the reverse direction to 23.2% — at or near chance. Masked-block
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
