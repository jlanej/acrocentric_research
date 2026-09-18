## 3. The compartments, one by one

### 3.1 Telomere and the distal p13 subtelomeric region

Each short arm terminates in a canonical telomere and a subtelomeric region containing the "ACRO" repeat family — a satellite shared across acrocentric distal ends that serves in practice as the marker confirming that an assembled contig really reaches the distal tip [@lin2026]. Distal p13 is otherwise a mosaic of satellite blocks and segmental duplications whose identity differs sharply between arms, and it is the single hardest compartment to assemble correctly: only 41% of short-arm scaffolds in the best-resourced pedigree study to date recovered both distal and proximal sequence, and only 20–40% of distal allelic sequence could be aligned between homologues in the same individual, with most aligned HSat1 showing only ~95% identity [@lin2026]. Distal sequence is where allelic divergence is highest, where assembly errors concentrate, and where the "same" region on two homologues may not be recognisably the same region at all.

This has an underappreciated consequence for variant analysis: **on the distal arm, the concept of an allele at a position is frequently undefined**, not because the sequence is hard to read but because there is no correspondence between the two haplotypes to define a position against. Section 13 develops this point.

### 3.2 The distal junction (DJ)

The DJ is the most important single feature of the acrocentric short arm for anyone trying to obtain a foothold in it. It is a ~300 kb segmental duplication, present in one copy on each of the five arms, immediately distal to the rDNA array and in inverted orientation relative to rDNA transcription [@floutsakou2013]. It is the most highly conserved sequence shared by all five arms; its discovery and characterisation defined the modern concept of the shared NOR architecture [@floutsakou2013,@mcstay2016].

Three properties make the DJ analytically central.

**It is transcribed and chromatin-marked.** The DJ carries promoter and active-transcription chromatin signatures that are conserved across cell types, and produces long non-coding transcripts [@floutsakou2013,@vansluis2019,@lin2026]. It is not inert heterochromatin: it is the functional interface between the rDNA array and the rest of the nucleus.

**It anchors the NOR to perinucleolar heterochromatin.** The DJ is positioned at the nucleolar periphery and is implicated in tethering the NOR within nucleolar architecture, in maintaining the boundary between transcribed rDNA and flanking heterochromatin, and in the reorganisation that follows rDNA damage [@mcstay2016,@vansluis2019]. Models of nucleolar assembly treat the DJ as the structural element that makes an rDNA array into a NOR.

**It is a usable anchor, and it varies.** Because the DJ is ~300 kb of non-satellite sequence with chromosome-distinguishing variants, it is the practical handle for identifying which arm a distal contig belongs to. Its total length is remarkably stable — 297–333 kb across a four-generation pedigree — yet a minigraph [@li2020minigraph] of twelve error-free DJ sequences resolved into twelve distinct bubbles, three of them nested, indicating extensive internal structural variation; a chr14 DJ differed from a chr22 DJ by a 1,525 bp deletion and three insertions [@lin2026]. Both facts matter: the DJ is conserved enough to be found reliably and variable enough to distinguish arms and haplotypes. **We regard the DJ as the natural coordinate origin for the distal arm** (section 22.1).

### 3.3 The rDNA array

The rDNA array is the reason the acrocentric short arms exist as a distinct biological category, and the reason they remain unresolved.

**Unit structure.** The human rDNA repeat unit is 44,838 bp in the reference sequence (GenBank KY962518), of which only 13,332 bp is the 47S pre-rRNA transcription unit: 5' external transcribed spacer, 18S rRNA (1,869 bp), internal transcribed spacer 1, 5.8S rRNA (157 bp), ITS2, 28S rRNA (5,051 bp), 3' ETS. The remaining 31.5 kb — 70% of the unit — is intergenic spacer carrying named tandem-repeat blocks (TR1, TR2), a long repeat (LR1), three simple-sequence-repeat blocks, the spacer promoter and enhancer region, and 26 Alu-family elements making up 12.4% of the unit [@ky962518]. Figure 4 maps this.

!FIG fig4_rdna.png | **Figure 4. The rDNA unit, the CHM13 rDNA models, and what the models are not.** (a) Feature map of the 44,838 bp human rDNA reference unit (GenBank KY962518). Only the leftmost 13.3 kb is transcribed by RNA polymerase I into 47S pre-rRNA; the intergenic spacer carries structured tandem-repeat and simple-sequence-repeat blocks and is dense in Alu elements. (b) Number of units in each CHM13v2.0 rDNA model, obtained by dividing modelled array length by the reference unit length. (c) These arrays are computational models; statistics computed inside them describe the model.

**Array organisation.** Units are arranged head-to-tail, uniformly oriented with transcription directed centromere-ward, in a single array per arm bounded distally by the DJ and proximally by the PJ [@nurk2022,@floutsakou2013,@mcstay2023]. Within an array, units are strikingly homogenised — the array behaves as a concerted-evolution unit, with sequence variants shared along a cluster rather than distributed at random [@hori2021,@agrawal2018]. Across arrays on different chromosomes, homogenisation is weaker but still present, which is the sequence signature of exchange between non-homologous arms.

**Copy number.** The five arrays together typically encode 200–600 rDNA genes per diploid genome, with wide interindividual variation [@nurk2022,@guarracino2023,@stults2008]. Copy number is coupled to expression and to mitochondrial abundance [@gibbons2014], varies in a concerted fashion so that total dosage is partially buffered across arrays [@gibbons2015], is unstable in cancer and with age [@stults2009,@xu2017,@wang2019rdna], and is heritable in a pattern that reflects array-level rather than base-level inheritance [@stults2008]. The CHM13 models contain 221 units haploid — within the expected range but not a measurement.

**Sequence variation.** rDNA units are not identical. Variant positions exist within the transcribed region, are conserved across populations, and are expressed tissue-specifically, implying functionally distinct ribosome populations [@parks2018]. The intergenic spacer is far more variable than the transcribed unit. Long-read work on chromosome 21 rDNA via TAR cloning first showed the extent of unit-to-unit variation directly [@kim2018], and *ribotin* now reconstructs distinct rDNA "morphs" — recurrent unit haplotypes — from HiFi and ONT data [@rautiainen2024ribotin]. Morphs, not units, are the natural allele for rDNA.

**Why it does not assemble.** Near-identical 45 kb units tiled across 0.7–3.6 Mb defeat every assembler: the repeat period vastly exceeds HiFi read length and approaches the practical limit of ultra-long ONT reads, so the graph contains a single collapsed loop with no way to count traversals. Even in the best current pipelines the rDNA array is a *predictable* collapse — an expected failure, annotated and excluded rather than resolved [@lin2026]. Our anchor analysis quantifies this: rDNA is 0.0% short-read accessible, 7.7% anchorable at 10 kb and 16.1% at 100 kb, with a median minimum-unique length of ~565 kb even within the artificially varied model.

### 3.4 The proximal junction (PJ)

The PJ flanks the rDNA array on the centromeric side, is also shared among arms, and is less well characterised than the DJ — in part because it sits adjacent to, and is partly composed of, satellite sequence, making it harder to isolate and harder to assemble [@floutsakou2013,@mcstay2023]. Functionally it is the second boundary element of the NOR, separating transcribed rDNA from proximal heterochromatin. For analysis it provides a second, weaker anchor into the arm; in practice proximal-arm contigs are more often anchored through the centromere and long arm than through the PJ.

### 3.5 β-satellite (Sau3A / BSR)

β-satellite is a GC-rich family built on a ~68 bp monomer, organised into higher-order arrays, and it is the satellite most specifically associated with acrocentric short arms — although it also appears at 1q and on the Y [@altemose2022sat]. It totals 6.12 Mb across the five short arms in CHM13 and is the dominant satellite of 22p (2.65 Mb, 21% of the arm). β-satellite is relatively tractable by satellite standards: 7.8% short-read accessible (the highest of any satellite compartment here apart from monomeric α-satellite) and 98.5% anchorable at 10 kb, reflecting substantial internal divergence between and within arrays.

### 3.6 HSat1A and HSat1B

HSat1 was historically "satellite 1" and is now split into two families with different composition and distribution [@altemose2022sat,@hoyt2022]. HSat1A is extremely AT-rich, built on a ~42 bp repeat, and forms the largest single satellite block on 13p: a 4.9 Mb array spanning the entire distal third of that arm, with 8.79 Mb across all five arms. HSat1B is smaller (1.09 Mb total) and concentrated on 15p and 22p.

Two observations matter analytically. First, HSat1A arrays are **hypervariable in length**: a single chr15 HSat1A array in one pedigree founder reached 1.64 Mb, six times its CHM13 length and 57% of that haplotype's distal sequence, with no assembly error [@lin2026]. Second, HSat1A is **mutationally hot** — it carried the largest share of de novo SNVs assigned to any single satellite family, including four in one paternal transmission on 13p, and a 6,763 bp de novo deletion on 22p [@lin2026]. The AT-rich composition is a plausible mechanistic contributor: the excess of A>C substitutions on the short arms localises to AT-rich repeat sequence, consistent with oxidative lesions such as 8-hydroxyadenine [@lin2026].

Despite their size, HSat1A arrays are *not* the analytical bottleneck: 98.8% of HSat1A positions have a unique anchor at 10 kb. Divergence within these arrays is high enough that long reads place them. Their difficulty is length variation, not ambiguity.

### 3.7 HSat2 and HSat3

HSat2 and HSat3 derive from related pentameric and related simple repeats and form the classical C-band heterochromatin of 1q, 9q, 16q and Yq12 [@altemose2022sat]. On the acrocentric short arms HSat3 is the single largest satellite compartment — 15.54 Mb, 23.5% of all short-arm sequence — and it is the defining feature of 15p, where an 8.0 Mb block occupies the proximal half of the arm (52.5% of that arm). HSat2 is essentially absent (2 kb).

HSat3 array length is highly variable between haplotypes: 1.69–6.78 Mb on chr15 within one family [@lin2026]. A diverged HSat3 block sitting between the main HSat3 array and HSat1B on 15p shares ~93% identity with the main array but is entirely distinct from adjacent HSat1B, and the two blocks differ in CpG methylation (~50% versus ~75%) — evidence that satellite blocks are epigenetically as well as structurally distinguishable [@lin2026]. HSat3 transcription is a documented feature of cellular stress and of cancer, and HSat2/3 RNA is implicated in tumour biology [@ting2011,@zhu2018,@bersani2015,@hoyt2022].

### 3.8 SST1 (NBL2) — the recombinogenic macrosatellite

SST1, also called NBL2, is a ~1.4 kb-unit macrosatellite present on 13p, 14p and 21p (and elsewhere in the genome). It has moved from obscurity to centrality in acrocentric biology over the past three years for two reasons.

**It marks the pseudo-homologous regions.** Pangenome-graph analysis of HPRC assemblies identified shared homology domains between non-homologous acrocentric arms, centred on SST1, whose positional-homology entropy is consistent with hotspots of ectopic recombination [@guarracino2023]. Long-read assembly of Robertsonian translocation carriers then identified a common breakpoint associated with the SST1 repeat, implicating these sequences directly in the formation of the commonest human constitutional rearrangement [@delima2025rob].

**It is epigenetically distinctive and clinically relevant.** SST1 arrays are consistently hypermethylated (~78% CpG) relative to flanking sequence (~48%) across 13p, 14p and 21p, and array lengths are bimodally distributed from 21 to 104 kb [@lin2026]. Somatic SST1 hypomethylation accompanies tetraploidisation in colorectal cancer [@gonzalez2021]. SST1 is therefore a locus where genotype (array length), epigenotype (methylation) and rearrangement risk converge — an unusually favourable target, discussed in section 22.6.

A caution on mechanism: the single ectopic chr13–chr21 recombination event captured directly in a pedigree did **not** break at SST1. It broke within a 630 kb segmental duplication of 99.53% identity located 1.6 Mb proximal to the SST1 array, at a cluster of PRDM9 binding motifs [@lin2026]. The honest current position is that SST1-proximal sequence is a recombination-prone *neighbourhood* containing multiple high-identity substrates, and that SST1 itself may be a marker of that neighbourhood rather than the exclusive substrate.

### 3.9 α-satellite: monomeric, divergent HOR, and active HOR

α-satellite is organised into higher-order repeats (HORs) built from ~171 bp monomers; an active HOR array recruits CENP-A and the kinetochore, while flanking divergent-HOR and monomeric α-satellite form the pericentromeric transition [@altemose2022cen,@altemose2022sat,@talbert2022]. Across the five short arms (as defined here, i.e. excluding the active array) α-satellite totals 5.54 Mb: 3.77 Mb monomeric, 1.70 Mb HOR and 0.06 Mb divergent HOR, the HOR component being the α-satellite islands within 15p and 22p.

The critical fact for cross-arm analysis is the **identity of the active arrays themselves**:

| Chromosome | Active HOR family | Active array in CHM13 |
| --- | --- | --- |
| 13p / 21p | S2C13/21H1L | 1.95 Mb / 0.33 Mb |
| 14p / 22p | S2C14/22H1L | 2.62 Mb / 2.92 Mb |
| 15p | S2C15H1L | 1.02 Mb |

Chromosomes 13 and 21 share a centromeric HOR family, and 14 and 22 share a different one. This was known from hybridisation data long before sequence was available [@choo1988], and it means the paralogy that confounds short-arm analysis extends through the centromere and does not stop at the arm boundary. Reads from a 13 centromere and a 21 centromere are mutually confusable; so are reads from 14 and 22. Centromere length and composition also vary substantially between haplotypes and populations [@logsdon2024cen,@miga2019,@sullivan2020], and chromosome 21 centromere size is asymmetric in families with Down syndrome [@mastrorosa2025]. Monomeric α-satellite is, at 17.6%, the most short-read-accessible compartment of the short arm — the one place where short reads retain partial traction.

### 3.10 Segmental duplications and the acrocentric "community"

The five short arms should be thought of as a single, partially shared sequence community rather than five independent loci. Figure 5 quantifies this.

!FIG fig5_crossarm.png | **Figure 5. Each short arm has a near-identical twin on another acrocentric chromosome.** (a) Percentage of the query arm (row) covered by a segmental duplication whose partner lies in the labelled arm (column), at ≥99% sequence identity; the diagonal reports duplications internal to the same arm. (b) Percentage of each arm covered by a duplication whose partner lies on a *different* acrocentric chromosome, at three identity thresholds. Calls are SEDEF segmental duplications on T2T-CHM13v2.0, which exclude satellite-dominated intervals; satellite arrays add further cross-arm identity beyond what is shown.

At ≥99% identity, 60.5% of 21p, 51.9% of 13p, 44.6% of 14p and 22.5% each of 15p and 22p lie in a duplication whose partner is a different acrocentric short arm. The pairwise structure recovers the classical cytogenetic pairing preferences: 21p–13p is the strongest relationship (53.6% of 21p has a ≥99% partner in 13p), followed by 21p–14p and 13p–14p; 15p is the most autonomous arm. Crucially, these are *non-satellite* duplications — the SEDEF call set excludes satellite-dominated intervals, and the median satellite fraction of the acrocentric–acrocentric pairs counted here is 0.26. The cross-arm identity contributed by the satellite arrays themselves sits on top of these figures.

This is the structural basis of three phenomena that recur throughout the rest of this report: natural hybrid short arms in the population [@guarracino2023]; Robertsonian translocations [@delima2025rob,@gerton2024]; and the intractability of read assignment. Segmental duplications genome-wide show elevated mutation and interlocus gene conversion [@vollger2023sdmut]; the acrocentric short arms are the extreme case of that phenomenon, and 14 of 103 short-arm de novo SNVs had an identical match elsewhere in the same repeat family, consistent with interlocus gene conversion rather than independent mutation [@lin2026].
