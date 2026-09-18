# Part II · Working with them
%SUB References, data generation, assembly, validation, alignment and annotation

## 7. References and coordinate systems

### 7.1 The choice of reference is the single largest determinant of what you can see

GRCh38 represents the acrocentric short arms as 66.85 Mb of N across the five cytogenetic p-arms, with 14.05 Mb of real sequence concentrated on 21p and 22p (Figure 2b). Reads originating from short-arm DNA in a GRCh38-aligned sample do not disappear; they mismap — preferentially into the fragments of short-arm sequence that *are* present, into the pericentromeric long arms, and into unplaced scaffolds. This generates spurious coverage peaks, spurious heterozygosity, and a well-documented class of false-positive variant calls that a complete reference removes [@aganezov2022].

!FIG fig2_reference_access.png | **Figure 2. Reference completeness, duplication content and callability.** (a) Composition of each short arm by satellite family (colours as Figure 1); the figure at the right of each bar is total arm length in Mb. (b) Sequence actually present in GRCh38 over the cytogenetic p-arm, against the complete T2T-CHM13v2.0 short arm. (c) Percentage of each arm inside a SEDEF segmental duplication, against the genome-wide figure of 7.3%. (d) Percentage of each compartment inside the T2T short-read accessibility mask, against the autosomal average of 86.9%; compartment sizes pooled across the five arms are given in the labels.

Switching to T2T-CHM13 recovers the sequence and, per the T2T analyses, improves variant calling genome-wide [@aganezov2022,@nurk2022]. It does not, however, make the short arms callable: our measurement is that 4.27% of short-arm sequence is inside the T2T short-read accessibility mask, against 86.9% of the autosomal genome. **A complete reference converts "invisible" into "visible but unreliable" — which is progress, but only if the unreliability is represented.**

### 7.2 Practical guidance on reference choice

- **Short-read cohorts already aligned to GRCh38.** Do not re-align an entire cohort solely for the short arms; the yield does not justify it. Instead, add the acrocentric regions to your excluded-region set explicitly (rather than relying on them being N), and treat any variant called within the 14.05 Mb of GRCh38 short-arm sequence as requiring orthogonal confirmation. Re-align if you also want the other T2T benefits — corrected false duplications, medically relevant genes, centromeres.
- **New short-read cohorts.** Align to T2T-CHM13v2.0 (with the HG002 Y for male samples, as distributed) and carry the accessibility mask as a first-class annotation. Report short-arm results only within the mask, and report the masked fraction.
- **Long-read cohorts.** Assemble. Aligning long reads to any single linear reference discards most of the information that made long reads worth generating in these regions; see section 9.
- **Graph references.** A pangenome graph is the correct representation in principle and is increasingly practical [@liao2023,@wang2022hpp,@garrison2024pggb,@hickey2023mc]. Current human pangenome graphs are, however, built from assemblies that are themselves incomplete across the short arms, so the graph inherits the gaps. Graph genotyping [@siren2021giraffe,@ebler2022pangenie] works well in the SD-rich non-satellite compartments and poorly inside long arrays. Section 22.2 proposes what an acrocentric-competent graph module would need.

### 7.3 Liftover and coordinate translation

Liftover between GRCh38 and CHM13 fails in the short arms by construction: there is no source interval to lift from. Even where GRCh38 has sequence, multi-mapping chains make the translation ambiguous. Two rules:

1. **Never lift short-arm coordinates.** Re-derive them by direct alignment or annotation transfer, and record which assembly and which annotation version produced them.
2. **Always record the annotation version.** CenSat v2.0 and v2.1 differ; RepeatMasker versions differ; the rDNA models are versioned. A short-arm coordinate without a reference *and* annotation version is not interpretable.

### 7.4 Masks you should be carrying

Four T2T-CHM13v2.0 resources should be attached to any short-arm analysis, and all four are used in this report:

| Resource | What it gives you | How we used it |
| --- | --- | --- |
| CenSat v2.1 satellite annotation | satellite family and array boundaries | composition, compartment stratification |
| SEDEF segmental duplications (`SD.full.bed`) | pairwise duplications with identity | cross-arm paralogy (Figure 5) |
| Short-read accessibility mask | regions where short-read genotypes are trusted | callability (Figures 1, 2d) |
| Minimum-unique-*k*-mer-length (MU) tracks | shortest genome-unique substring at each position | anchor ladder, regime map (Figures 3, 6) |

The MU tracks in particular are underused. They answer, per base, "how much contiguous sequence do I need to place a read here uniquely?" — which is precisely the question that determines technology choice, and is more informative than a single-read-length mappability track.

## 8. Data generation: what to buy for what question

Figure 3 is the basis for the recommendations in this section.

!FIG fig3_anchor_ladder.png | **Figure 3. Unique-anchor availability as a function of read length, and where the hard sequence is.** (a) Percentage of positions in each short arm at which the shortest genome-unique substring is no longer than the abscissa — equivalently, the fraction of positions at which a read of that length can be placed uniquely in T2T-CHM13v2.0. The grey reference curve is a 20 Mb euchromatic interval of chr13q. Shaded bands mark the working read lengths of current platforms. (b) Median minimum-unique-anchor length in 10 kb windows along each arm; red bars mark the rDNA arrays.

### 8.1 The read-length ladder, read as a purchasing decision

At 150 bp, 15.6% (13p) to 45.4% (14p) of short-arm positions are uniquely anchorable, against 98.9% in euchromatin. At 1 kb — the scale of linked reads or short-insert mate pairs — this rises to 49.3–84.7%. At 10 kb, the HiFi working range, 77.5–93.6%. Then the curve flattens: from 10 kb to 100 kb, gains are 1.5–3.4 percentage points. From 100 kb to 1 Mb, a second rise appears on 13p (79.4%  to  86.1%) and 21p (80.1%  to  94.1%) as the rDNA models become traversable.

The purchasing implications are unusually clear:

- **Short reads alone**: adequate for the monomeric α-satellite and parts of β-satellite and the non-satellite compartment, within the accessibility mask; useless for arrays. Do not buy more short-read depth expecting to reach the arms — the limit is uniqueness, not coverage.
- **HiFi (~15–25 kb)**: the single highest-value purchase. It captures essentially the whole non-rDNA arm. Depth matters: published short-arm assembly used ~40× HiFi [@lin2026].
- **Ultra-long ONT (>100 kb)**: buys little additional *anchoring* but is essential for *assembly* — it spans array-internal repeat periods and resolves graph tangles that HiFi cannot. ~30× UL-ONT with N50 above 100 kb is the working figure [@lin2026,@antipov2025verkko2].
- **Hi-C / Micro-C (~30×)**: the decisive addition for the distal arm. Proximity ligation is what scaffolds distal to proximal sequence across the unassemblable rDNA array; it raised distal+proximal recovery to 41% of short arms in the best current study [@lin2026,@antipov2025verkko2]. Without it, distal sequence is usually orphaned.
- **Strand-seq**: independent structural validation and haplotype phasing, particularly for inversions and misorientations [@sanders2020,@porubsky2022inv].
- **Optical mapping**: still the most direct way to measure multi-megabase array length, and a useful orthogonal check on assembly-derived array sizes.
- **Targeted approaches**: Cas9-guided nanopore enrichment [@gilpatrick2020] makes single-locus questions — one rDNA array, one SST1 array — affordable at cohort scale, and is badly underused in this region.
- **Single-molecule chromatin state**: chromatin fibre sequencing resolves regulatory architecture per molecule [@stergachis2020] and is one of the few ways to ask whether individual rDNA units or satellite blocks are active in a given cell. **Methylation**: obtained free with ONT or HiFi and genuinely informative here, since CDR position, SST1 hypermethylation, satellite block boundaries and rDNA activity are all methylation-visible [@gershman2022,@lin2026].

:::note A DNA-extraction caveat that invalidates results
Ultra-long protocols require ultra-high-molecular-weight DNA, and the acrocentric short arms are exactly where fragment length is limiting. Phenol-chloroform or dedicated HMW kits from fresh cells are standard [@lin2026]. Equally important: lymphoblastoid cell lines accumulate somatic changes, and satellite arrays and rDNA copy number are among the first things to drift. For mutation or copy-number work, primary tissue is strongly preferable — the pedigree study excluded two generations specifically because only cell lines were available [@lin2026].
:::

### 8.2 Minimum viable designs

| Question | Minimum data | Expected outcome |
| --- | --- | --- |
| Is this individual a Robertsonian carrier? | 20× HiFi, or karyotype | whole-arm fusion detectable from assembly graph topology or contig structure |
| rDNA copy number per individual | 30× short-read WGS + a calibrated estimator, or ddPCR | total copy number ±10–15%; not per-chromosome |
| rDNA morph composition | 30× HiFi (+ ONT) and `ribotin` | morph catalogue and relative abundance; not array order |
| Satellite array lengths, one arm | 40× HiFi + 30× UL-ONT | array length to ~kb where the array assembles |
| Complete phased short arms | 40× HiFi + 30× UL-ONT + 30× Hi-C | ~90% of expected arm length, rDNA collapsed, ~41% with distal scaffolded |
| De novo SNVs on short arms | the above, for both parents and child | ~30 Mb callable per haplotype; 1–10 SNVs per transmission |
| Cohort-scale array-length association | short reads + *k*-mer or depth estimators | relative array-size estimates, cohort-normalised |

## 9. Assembly

### 9.1 What current assemblers actually achieve

The path to this point ran through single-chromosome proofs of principle [@miga2020x], the first complete haploid genome [@nurk2022], and semi-automated diploid pipelines [@jarvis2022]; the general capabilities and limits of long-read human sequencing are reviewed elsewhere [@logsdon2020].

The honest summary is that **complete short arms are not yet routine, and rDNA arrays do not assemble at all.** Typical diploid long-read assemblies retain 50–100 gaps, and many map to the acrocentric short arms [@porubsky2023gaps]. With 40× HiFi, 30× UL-ONT and 30× Hi-C and a purpose-built assembler, the current state of the art is: ~90% of expected short-arm length recovered per haplotype; rDNA arrays predictably collapsed; 4.8% of bases misassembled or collapsed (~0.72 Mb per haplotype per chromosome); distal and proximal sequence successfully scaffolded together in 41% of arms; 97% of expected haplotypes recovered as pq-scatigs when rDNA and distal sequence are not required [@lin2026].

Plan around these numbers rather than against them. An analysis design that requires complete, error-free short arms will fail; one that requires 30 Mb of *validated callable* sequence per haplotype will succeed.

### 9.2 Assembler choice

**Verkko2** is the current tool of choice for this problem. It builds long-read de Bruijn graphs from HiFi plus ONT and integrates proximity-ligation data for phasing and scaffolding — capabilities added specifically to improve acrocentric short-arm assembly and to handle rare translocations including Robertsonians [@antipov2025verkko2,@rautiainen2023verkko]. **hifiasm** in UL mode is the main alternative and is faster [@cheng2021hifiasm,@cheng2024hifiasmul]; for centromeres and short arms, published short-arm work has used Verkko2 [@lin2026]. Polishing with DeepPolisher improves base accuracy [@mastoras2025]; expect QV in the mid-50s [@lin2026].

For the rDNA specifically, **`ribotin`** is the appropriate tool: rather than attempting a linear assembly, it reconstructs and phases rDNA morphs from reads that map to the rDNA consensus, producing a morph catalogue and abundance estimates [@rautiainen2024ribotin]. This is the right shape of answer for the rDNA, and section 22.3 argues it should be generalised to other arrays.

### 9.3 Identifying and anchoring short-arm contigs

A short-arm contig is worthless unless you can say which chromosome it came from. The established procedure anchors through the long arm [@guarracino2023,@lin2026]:

1. Align assembly contigs to T2T-CHM13 with `minimap2 -x asm20 --secondary=no -s 25000 -c --eqx --cs` [@li2018mm2].
2. Discard alignments shorter than 100 kb or below 90% identity.
3. Require ≥1 Mb of q-arm alignment, at least 1 Mb away from the centromere, so that the anchor is in sequence that segmental duplications do not make ambiguous.
4. Retain the contig as a pq-scatig spanning the p-arm, centromere and 5 Mb of q-arm.
5. Confirm distal completeness by the presence of the subtelomeric ACRO repeat and a DJ flanking the rDNA array.

Steps 3 and 5 are where most errors enter. A p-arm sequence anchored only through a short or centromere-proximal q-arm alignment can be assigned to the wrong chromosome, and 13/21 and 14/22 are the pairs on which this happens.

## 10. Validation and quality control

Short-arm assemblies fail in characteristic ways, and every one of them produces plausible-looking variant calls. The following are not optional.

**Collapse and misassembly detection.** `NucFreq`/`NucFlag` flags collapsed repeats and misjoins from the frequency of the second-most-common base in read pileups against the assembly — collapse shows up as a region where the second base count exceeds a threshold [@vollger2022sd]. `HMM-Flagger` uses a hidden Markov model on read coverage to classify erroneous, falsely duplicated and collapsed intervals [@liao2023]. Use both and take the union as the excluded set; they have different failure modes [@lin2026].

**Base-level accuracy.** `Merqury` compares assembly *k*-mers against short-read *k*-mers to estimate QV and phasing accuracy without a reference [@rhie2020merqury].

**Structural visualisation.** `ModDotPlot` for self-similarity dotplots of each arm — the fastest way to see a collapsed array or an inverted misjoin [@sweeten2024moddotplot]. `SVbyEye` for all-versus-all comparisons across haplotypes and family members [@porubsky2025svbyeye].

**Centromere-specific QC.** `CenMAP` identifies centromeric contigs, runs NucFlag, and annotates HOR organisation with RepeatMasker and HumAS-HMMER [@lin2026,@altemose2022cen].

**Orthogonal read support.** For any variant, require support from two independent long-read technologies in the child and absence in both parental read sets [@lin2026]. This is the single most effective false-positive control available in these regions, and it is the reason the pedigree DNM estimate is credible.

**General polishing and validation strategy.** The T2T consortium's validation and polishing framework remains the reference treatment [@mccartney2022].

:::note The QC step people skip, and why it matters most
Because similar arrays occur on both the distal and proximal sides of the rDNA within the same arm, an assembly that is missing distal sequence will happily align its *proximal* array to the *distal* array of a relative — producing a confident, entirely false, apparent difference. The fix is a 10 kb sliding-window identity check across every parent–offspring alignment, excluding any window below 99.9% identity [@lin2026]. Without it, distal/proximal misalignment is a systematic and invisible source of false de novo calls.
:::

## 11. Alignment, and the assignability problem

### 11.1 Mappability is not assignability

Every mappability metric, including the MU tracks used here, answers a question about the *reference*: is this substring unique in this assembly? The question that matters for analysis is different: given a read from *this individual*, can it be assigned to the correct locus and the correct haplotype? Call the first *mappability* and the second *assignability*.

In euchromatin the two coincide. In the acrocentric short arms they diverge sharply, for three reasons:

1. **Length variation destroys positional correspondence.** If a sample's HSat1A array is 1.64 Mb and the reference's is 0.27 Mb [@lin2026], reads from the sample's array will map somewhere within the reference array with high confidence and high identity, and the position will be meaningless.
2. **Paralogy is individual-specific.** Which arm a given block is most similar to depends on the individual's haplotypes, including whether they carry a hybrid arm [@guarracino2023]. A read that is uniquely placeable against CHM13 may be genuinely ambiguous against the sample's own genome.
3. **The reference is one haplotype of a hypervariable region.** CHM13 is a single hydatidiform-mole-derived haplotype. For the short arms it is not a consensus in any meaningful sense.

The practical consequence: **MAPQ is not a measure of correctness here.** High MAPQ against CHM13 in a satellite array is frequently high confidence in a wrong answer. Conversely, MAPQ 0 correctly reports ambiguity and should not be filtered away silently — it should be counted and reported.

### 11.2 Strategies that work

**Map to a personalised reference.** The highest-value single change. In trio work, align the child's reads to the *parental haplotype assemblies* rather than to a population reference; de novo variants then appear as differences against the actual template from which the child's chromosome was copied [@lin2026]. This converts an intractable multi-mapping problem into a tractable two-haplotype problem, and it is what made short-arm DNM measurement possible. The same logic applies outside trios: align to the sample's own assembly.

**Use paralog-specific variants (PSVs).** Where two paralogous arms must be distinguished, the informative signal is the set of positions at which they differ. Extracting PSVs from a multiple alignment of the candidate paralogs, then genotyping reads at those positions, is the established approach for SD-rich regions and is how the ectopic recombination breakpoint was refined to 18.7 kb [@lin2026,@porubsky2025_22q].

**Align assembly to assembly, not reads to reference.** All-versus-all alignment of haplotype assemblies with `minimap2 -x asm20 -c --eqx -D -P --dual=no`, or `wfmash` for larger sets, followed by sliding-window identity computation, is the workhorse comparison in this region [@li2018mm2,@garrison2024pggb,@lin2026].

**Use graphs where they fit.** Pangenome graphs handle the SD-rich non-satellite compartments well and make paralogy explicit rather than implicit [@garrison2024pggb,@hickey2023mc,@guarracino2022odgi]. They currently handle megabase tandem arrays poorly.

**Report the denominator.** Whatever strategy is used, the callable denominator must be computed and reported per sample and per arm. An unstated denominator is the commonest defect in short-arm claims.

## 12. Annotation

Annotation of a new short-arm assembly is a distinct step from assembly and should be done explicitly.

**Satellite annotation.** The current practical approach maps a curated set of target satellite consensus sequences — derived from CenSat v2.1 plus the HSat consensus set — onto the assembly, inverting the usual query/reference orientation so that all matches to each target are recovered rather than only the best one [@lin2026,@altemose2022sat,@censat]. Alignment blocks within 500 bp are merged and filtered by length to suppress spurious LINE matches. This is a pragmatic procedure, not a principled classifier, and the field would benefit from a proper probabilistic satellite annotator (section 20.9).

**α-satellite HOR structure.** `HumAS-HMMER`/HumAS-SD for HOR monomer classification, integrated in `CenMAP` [@altemose2022cen,@lin2026].

**Repeats and masking.** RepeatMasker, Tandem Repeats Finder and WindowMasker in combination, with the union used as the soft-mask, before SD detection [@benson1999trf,@morgulis2006,@lin2026].

**Segmental duplications.** `SEDEF` on the repeat-soft-masked assembly, retaining duplications above 1 kb and 90% identity; the conventional additional filter excludes duplications more than 70% satellite [@numanagic2018sedef,@lin2026]. Be aware that this filter is what makes SD statistics and satellite statistics non-overlapping; it also means published SD content *understates* total paralogy.

**Landmarks.** Telomeres and gaps with `seqtk telo`/`gap`; DJ and PJ by alignment to reference junction sequences; rDNA by alignment to the KY962518 unit, rotated to begin upstream of the 45S transcription start site so that the promoter is included [@lin2026,@ky962518].

**Methylation.** Basecall with modification tags, align to the sample's own assembly, and pile up with `modkit`; 10 kb bins are a reasonable default for array-scale comparison [@lin2026].
