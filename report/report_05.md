# Part III · Variants, from single nucleotides to whole arms
%SUB A framework, then class-by-class practice

## 13. A three-regime framework, and why VCF breaks

### 13.1 The framework

Compartment labels (HSat3, β-satellite, rDNA) describe sequence content. They are not the right axis for deciding how to analyse a region, because compartments with the same name behave differently depending on array length and internal divergence — HSat1A is a satellite and is 98.8% anchorable at 10 kb; rDNA is a satellite and is 7.7% anchorable at 10 kb. The right axis is **how much contiguous sequence is needed to place a read uniquely**, which is exactly what the MU tracks measure. Figure 6 partitions the short arms on that axis into three regimes.

!FIG fig6_regimes.png | **Figure 6. Three analysis regimes across the acrocentric short arms.** (a) Every 10 kb window of each arm classified by its median minimum-unique-anchor length. (b) The proportion of each arm falling in each regime. (c) What each regime permits analytically. Regime boundaries at 1 kb and 100 kb correspond approximately to the transition from short-read to long-read addressability, and from long-read to array-scale-only.

**Regime I — read-addressable (unique anchor ≤ 1 kb).** 49–85% of each arm. Per-base genotypes are meaningful. Short reads recover part of this regime (the part with anchors under ~150 bp, plus paired-end reach); HiFi recovers essentially all of it. Standard SNV and small-indel calling applies, subject to the caveats in section 14.

**Regime II — long-read addressable (1–100 kb).** 10–40% of each arm; 39.7% of 15p. Per-base genotypes are obtainable, but only from assemblies or from reads that span the local repeat period, and the coordinates require haplotype context to be interpretable. Standard read-to-reference pipelines produce confident wrong answers here.

**Regime III — array-scale only (> 100 kb).** 4.8–20.6% of each arm; 20.6% of 13p and 19.9% of 21p. No per-base coordinate in this regime is portable between individuals. The meaningful observables are array-level: length, copy number, unit composition, methylation state.

### 13.2 Why the VCF model fails in Regimes II and III

A VCF record asserts that at position *p* of reference *R*, this sample carries allele *A* instead of the reference allele. That assertion presupposes that position *p* has a counterpart in the sample's genome. Inside a satellite array whose length differs between reference and sample by hundreds of kilobases — the normal case — there is no such counterpart. The aligner will nonetheless produce a position, and the caller will nonetheless produce a record, and the record will be reproducible, and it will still be meaningless: it describes the alignment's arbitrary choice of where to place the indel burden, not a property of the sample.

This is not a pedantic objection. It has three concrete consequences that we have seen repeatedly in practice and that are visible in the numbers above:

1. **Non-comparability across samples.** Two samples with identical arrays but different alignment paths get different VCF records. Two samples with genuinely different arrays may get the same record.
2. **Non-comparability across references.** The same sample re-aligned to a different assembly, or to its own assembly, yields a different record set — with no mapping between them.
3. **Accumulating false discovery in cohort analysis.** Because the errors are reference-driven, they are *correlated across samples* in ways that mimic real population structure. Any association test run over Regime II/III VCF records is testing alignment behaviour.

Our position is that **per-base variant representation should be abandoned in Regime III and used only with haplotype context in Regime II.** Section 22.4 sketches what should replace it.

### 13.3 The decision matrix

Table 6 maps variant class against regime to a recommended method. It is the operational summary of this report.

| Variant class | Regime I | Regime II | Regime III |
| --- | --- | --- | --- |
| **SNV / small indel** | standard calling; require accessibility mask for short reads; HiFi assembly-based calls preferred | assembly-to-assembly calling against a personalised reference; orthogonal two-technology read support mandatory | not attempted; report as not assessed |
| **Indel 50 bp – 10 kb** | long-read SV calling (`sniffles2`, `delly`), assembly confirmation | assembly-to-assembly with `Truvari` matching at 80% size/sequence; manual review | not attempted |
| **Tandem-repeat unit change (e.g. HOR *n*-mer gain/loss)** | tandem-repeat genotyping (`TRGT`) where the array is short and anchored | assembly-based unit-structure comparison; HOR annotation diff | array-level unit-spectrum comparison only |
| **Satellite array length** | direct from assembly | direct from assembly where the array assembles; otherwise optical mapping or read-depth | read-depth or *k*-mer estimate; optical mapping; no coordinates |
| **rDNA copy number** | — | — | calibrated read-depth, *k*-mer, or ddPCR; total and not per-chromosome |
| **rDNA sequence variation** | — | — | morph reconstruction (`ribotin`); report morphs and abundances, not positions |
| **Large SV / NAHR (> 100 kb)** | assembly all-vs-all; PSV-based breakpoint refinement | same, plus Strand-seq and Hi-C confirmation | inferred from flanking Regime I/II anchors; array content not resolvable |
| **Whole-arm exchange (Robertsonian, hybrid arm)** | assembly graph topology; contig-level pq structure; karyotype remains valid | same | same |
| **Methylation / CDR state** | standard modified-basecall pileup | same, on personalised reference | array-level mean methylation; CDR position relative to array start |

!PAGEBREAK

## 14. SNVs and indels, including de novo calling in trios

This section is written for the case we expect to be most common: a trio or pedigree cohort in which the short arms have been excluded, and the question is what can be recovered.

### 14.1 What is callable, and how to define it

Three nested denominators must be distinguished, and reported:

- **Reference sequence present**: 66.07 Mb (T2T-CHM13), 14.05 Mb (GRCh38).
- **Short-read accessible**: 2.82 Mb — 4.27% of the short arms. This is the honest short-read denominator.
- **Long-read/assembly callable**: ~30.7 Mb per haplotype in the best current pedigree analysis, after excluding rDNA, assembly-error regions flagged by NucFreq and Flagger, misaligned windows below 99.9% identity, and unanchored distal sequence [@lin2026].

The order-of-magnitude difference between the second and third figures is the entire argument for long reads in this region.

The callable-region definition that we recommend adopting, following [@lin2026], is the intersection of: (i) transmitted/aligned blocks ≥1 Mb at ≥99% identity to the template haplotype; (ii) outside the union of NucFreq and Flagger flags on *both* query and reference; (iii) outside any 10 kb window below 99.9% identity between parent and child; (iv) outside rDNA models and unresolved arrays. Report it per sample, per haplotype, per chromosome.

### 14.2 De novo SNV calling: the procedure that works

The assembly-first, personalised-reference procedure is currently the only one with demonstrated validity in these regions [@lin2026]:

1. Assemble and phase both parents and the child (40× HiFi, 30× UL-ONT, 30× Hi-C each); anchor short arms as pq-scatigs.
2. Align child haplotypes to parental haplotypes (`wfmash -s 50k -l 150k -p 90 -n 1`), keeping blocks >1 Mb at >99% identity as the transmitted set.
3. Call variants from the haplotype-to-haplotype alignment (`wgatools call`) [@wei2025wgatools].
4. Require every candidate to be outside assembly-error regions on both sequences.
5. Require ≥3 supporting child reads in *both* HiFi and UL-ONT, and absence from both parental HiFi and UL-ONT read sets.
6. Manually review every candidate.
7. Where a further generation exists, use transmission as validation — 85% of G3 candidates were confirmed as transmitted [@lin2026].

Expect 1–10 de novo SNVs per transmission per individual across all five arms; for comparison, familial long-read sequencing raises de novo yield genome-wide relative to short-read trio calling [@noyes2022]; ~80% will fall in annotated satellite or SD sequence [@lin2026].

### 14.3 What fails, and the specific false-positive classes

- **Short-read trio DNM calling in the short arms is not salvageable.** With 4.27% accessibility and pervasive paralogy, the false-positive rate from mismapping exceeds any plausible true rate. If you have only short reads, the correct action is to exclude the short arms *explicitly and with a stated denominator*, not to filter aggressively and hope.
- **Distal/proximal misalignment** produces false de novo calls with full read support on both sides, because the reads genuinely support the sequence — just at the wrong locus. Only the sliding-window identity check catches it.
- **Collapsed arrays** produce apparent heterozygosity and apparent de novo variants at high read support. NucFreq is the control.
- **Interlocus gene conversion masquerades as de novo mutation.** Fourteen of 103 short-arm de novo SNVs matched an existing allele elsewhere in the same repeat family [@lin2026]. These are real sequence changes but the mechanism is conversion, not point mutation, and they should be annotated separately: check every candidate against the set of alleles present elsewhere in the same repeat family in the same individual.
- **Cell-line artefacts.** Satellite and rDNA changes accumulate in culture. Use primary material for mutation work.

### 14.4 Recalibrating expectations

If you do recover short-arm DNMs, the expected rate is ~10× the autosomal rate, with depleted CpG>TpG and elevated C>G and A>C [@lin2026]. Two practical consequences: a short-arm DNM count that looks like the autosomal rate is probably an underestimate reflecting incomplete ascertainment, not a negative result; and a mutational-signature analysis that pools short-arm and euchromatic DNMs will be distorted by the short-arm spectrum.

## 15. Tandem-repeat and array-level variation

Array-level variation is the dominant mode of variation in these regions, and it is the class least well served by existing tooling.

**What to measure.** Array length; unit copy number; unit-type composition (for HOR arrays, the *n*-mer spectrum); internal divergence; boundary positions relative to flanking anchors; methylation profile.

**From assemblies.** Where an array assembles cleanly, measure it directly and report the flanking anchor sequences used to delimit it. Validate with NucFreq — a collapsed array reports a spuriously short length, which is the commonest error in published array sizes.

**Unit-level SVs.** De novo unit gains and losses in α-satellite HOR arrays are detectable and real: stepwise deletion of a 1,364 bp 8-mer on chr22 in two independent transmissions, and insertion of a 680 bp 4-mer on chr21 [@lin2026]. These are best represented as changes in the HOR unit spectrum rather than as indels at coordinates. `TRGT` is the appropriate genotyper where arrays are short and anchorable [@dolzhenko2024trgt]; the broader population context for tandem-repeat variation and its expression consequences provides the comparison baseline [@trvar].

**From short reads.** Relative array length can be estimated from read depth over array-specific *k*-mers, with cohort-internal normalisation. This is a *relative* measure, useful for cohort comparisons and association testing, and it does not yield per-chromosome resolution for families shared between arms. Validate against a subset with long-read assemblies.

**Report array length as a quantitative trait**, with units, the estimator, and the normalisation. This is the representation that survives comparison across studies; a VCF indel record is not.

## 16. rDNA copy number and sequence

The rDNA deserves separate treatment because it is simultaneously the most biologically consequential and least tractable compartment.

**Copy number.** Diploid totals are typically 200–600 units with wide variation [@nurk2022,@stults2008]. Estimation options, in decreasing order of practicality:

- *Short-read depth*: coverage over the rDNA unit relative to genome-wide coverage, with GC and mappability correction. Cheap, cohort-scalable, sensitive to library preparation and to the reference unit used. Requires within-cohort normalisation and, ideally, calibration against an orthogonal method on a subset. Note that the CHM13 rDNA models must be masked; use the KY962518 unit as the target instead.
- *ddPCR*: the accuracy standard for total copy number; low throughput.
- *k-mer based*: count rDNA-specific *k*-mers in raw reads; avoids alignment entirely and is robust, and this is the approach we would recommend for large cohorts.
- *Long-read depth over collapsed arrays*: usable, but long-read depth is less uniform.
- *Optical mapping / FISH*: per-chromosome array size, low throughput, and the only routes to chromosome-resolved counts.

**Per-chromosome assignment is not currently solvable at scale**, and claims of per-chromosome rDNA copy number from short reads should be treated sceptically. Chromosome-specific IGS variants exist and are the plausible basis for a solution (section 22.5).

**Sequence variation and morphs.** Use `ribotin` to reconstruct morphs from HiFi (+ONT) [@rautiainen2024ribotin]. Report the morph catalogue and abundances. Transcribed-region variants are functionally interesting and tissue-specifically expressed [@parks2018]; IGS variants are more numerous and are the better barcodes. Do not report rDNA variants as VCF records against a single consensus unit: the same variant in a different unit is a different allele, and abundance rather than genotype is the meaningful quantity.

**Methylation and activity.** The active/silent unit ratio is regulated and is not inferable from copy number [@grummt2003]. Long-read methylation over rDNA provides a proxy and should be reported alongside copy number in any study proposing rDNA dosage as a phenotype.

## 17. Large structural variants, NAHR and whole-arm events

**Ectopic exchange and hybrid arms.** Detection requires assemblies and all-versus-all comparison; the diagnostic signature is a contig whose alignment transitions from one chromosome's haplotype to another's at a defined point, with read coverage from the sample dropping at the breakpoint on one template and resuming on the other [@lin2026]. Breakpoint refinement uses paralog-specific variants in a multiple alignment of the candidate paralogs, reaching ~18 kb resolution in the published case. Confirm with Strand-seq and, where available, Hi-C contact patterns.

**Robertsonian translocations.** Complete sequencing is now possible and has been done [@delima2025rob]; Verkko2 explicitly supports translocation assembly [@antipov2025verkko2]. For detection at scale, the assembly graph is more informative than any alignment: a whole-arm fusion appears as a contig carrying two different long arms joined through acrocentric short-arm/centromeric sequence. Karyotype remains a valid and cheap orthogonal test and should not be abandoned.

**Inversions.** Pericentromeric inversion polymorphisms are common on acrocentric long arms — a 515 kb 14q11 inversion segregates in the CEPH1463 pedigree, and q-arm inversion polymorphisms were found on chr14, chr15 and chr22 [@lin2026,@porubsky2022inv]. They matter for short-arm work because they delimit recombination and can disrupt synapsis; the three chr15 crossovers observed in that pedigree all mapped immediately adjacent to an SD-mediated inversion [@lin2026]. Strand-seq is the method of choice [@sanders2020].

**Isodicentrics and marker chromosomes.** Acrocentric-derived small supernumerary marker chromosomes are a recurrent clinical finding, usually characterised cytogenetically. Sequence-level characterisation is now feasible and largely undone.

## 18. Somatic and cancer applications

Three observables are within reach and largely unexploited.

**rDNA copy-number change in tumours.** Measurable from existing short-read tumour/normal WGS with a calibrated estimator, at no additional sequencing cost. Given the established instability of rDNA in cancer [@stults2009,@xu2017,@wang2019rdna], this is an obvious retrospective analysis on existing cohorts.

**Satellite transcription and hypomethylation.** HSat2/3 transcription and SST1 hypomethylation are established cancer-associated changes [@ting2011,@bersani2015,@zhu2018,@gonzalez2021]. Both are measurable from data types that large tumour cohorts already have, provided a complete reference is used.

**Mosaic and clonal structural change.** Acrocentric whole-arm events and satellite array changes in somatic tissues are essentially uncharacterised. Sensitive mosaic SV detection is available for long reads [@smolka2024], and read-based callers and integrated pipelines developed for short and long reads [@rausch2012delly,@kolmogorov2023,@behera2024dragen] can be applied within the anchorable fraction and, applied to acrocentrics, would answer whether somatic acrocentric instability is a real phenomenon at appreciable frequency.

## 19. Short-read cohorts: what is actually possible

Most existing cohorts — including most rare-disease and pediatric cohorts — are short-read. The following is what we consider defensible, and what we consider not.

**Defensible.**

1. **Explicit exclusion with a stated denominator.** Define the acrocentric short arms in your reference's coordinates, exclude them from variant analysis, and state the excluded base count in methods. This is strictly better than relying on N-masking, which silently redistributes the reads.
2. **Accessibility-masked calling on T2T.** Within the 2.82 Mb accessible fraction, short-read genotypes are usable. Report results only there, with the mask cited.
3. **Total rDNA copy number as a quantitative trait.** Cohort-normalised, *k*-mer or depth-based, validated on a long-read subset. This is a real, heritable, biologically consequential phenotype [@gibbons2014,@gibbons2015] and it is obtainable from data you already have.
4. **Relative satellite array-size estimates.** Family-specific *k*-mer depth for HSat1A, HSat3, β-satellite and SST1, cohort-normalised. Useful as quantitative traits and for outlier detection.
5. **Robertsonian carrier screening by an unusual-coverage signature.** A ROB carrier has lost most of two short arms; the resulting deficit in short-arm satellite and rDNA *k*-mer depth is a candidate screening signal. This needs validation against known carriers (section 22.6), but the data requirement is nil.
6. **Outlier detection rather than genotyping.** Even where per-base genotypes are meaningless, gross per-sample summaries (total satellite *k*-mer content per family, rDNA copy number, arm-level coverage profile) are comparable within a cohort and will flag the extreme individuals who merit long-read follow-up. We regard this as the highest-value short-read strategy available.

**Not defensible.**

- Per-base SNV or indel calls outside the accessibility mask.
- De novo mutation calls anywhere in the short arms from short reads alone.
- Per-chromosome rDNA copy number.
- Structural-variant calls inside satellite arrays.
- Any claim of *absence* of variation in these regions. With 4.27% accessibility, absence of calls carries essentially no information — and this is the most frequent unforced error in the literature. The correct phrasing is "not assessed", not "no variants found".
