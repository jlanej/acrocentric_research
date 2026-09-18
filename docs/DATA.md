# Data provenance

Every input is public and small; no genome assembly is downloaded.
`analysis/fetch_data.sh <work>` retrieves all of it (~105 MB).

## Per-haplotype acrocentric annotation

| File | Source |
| --- | --- |
| `all_samples_censat.bed` | Platinum Pedigree Consortium, `AcroMutRecomb` repository, `annotation/`. Satellite blocks for 161 acrocentric arms of the CEPH1463 pedigree. Carries ACRO, DJ, PJ and SST1 labels that the standard CenSat vocabulary does not. |
| `flagger_nucfreq_merged.bed` | same repository. Union of Flagger- and NucFreq-flagged intervals; the assembly-error exclusion set for those arms. |
| `hprc/<sample>_hap<N>_hprc_r2_v1.cenSat.bed` | HPRC release 2, `s3://human-pangenomics/working/HPRC/<sample>/assemblies/release2/annotation/censat/`. 410 haplotypes from 205 individuals. |
| `aux/<sample>_hap<N>_hprc_r2_v1.0.1.fa.gz.fai` | same prefix; scaffold lengths, used for chromosome assignment. |
| `aux/<sample>_hap<N>_hprc_r2_v1.active.centromeres.bed` | same prefix; active centromere intervals. |

Of the 410 release 2 haplotypes, 188 (94 individuals) yielded at least one
chromosome-scale acrocentric scaffold and enter the analysis; the rest are not
assembled to chromosome scale across the acrocentrics.

## T2T-CHM13v2.0 reference annotation (report figures 1-6)

All from `s3://human-pangenomics/T2T/CHM13/assemblies/annotation/`:
`chm13v2.0_censat_v2.1.bed` (satellite annotation), `chm13v2.0_SD.full.bed` (SEDEF
segmental duplications with identity), `accessibility/combined_mask.bed.gz`
(short-read accessibility mask) and `reference_accessibility_comparison.txt`,
`mappability/chm13v2.0.mul.bw` (minimum unique k-mer length; 3.2 GB, not fetched by
default), `chm13v1.1.rdna_model.bed` (rDNA model intervals - mask these),
`chm13v2.0_cytobands_allchrs.bed`. GRCh38 `gap` and `cytoBand` tables from UCSC.
The rDNA unit reference is GenBank `KY962518.1`.

## Published values used but not recomputed

De novo mutation counts and rates, recombination counts, DJ and SST1 length ranges,
and assembly-completeness figures are from Lin et al., *Cell* 189:4876-4890.e8 (2026),
doi:10.1016/j.cell.2026.05.035. Denominators for the per-compartment rates in report
section 5.4 were computed here from the annotation above.
