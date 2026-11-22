# Changelog

## 2.0.0 — 2026-11-22

Complete rewrite. The 1.x files were never loaded by Claude Code (flat `.claude/skills/*.md` plus an invented `settings.json` `skills` key); this release moves to the supported `.claude/skills/<name>/SKILL.md` layout with YAML frontmatter.

### Structure
- One directory per skill, `SKILL.md` ≤ 400 lines, long material in supporting files (`reference.md`, `bulk.md`, `atac.md`, …).
- Shared operating rules in the `genomics` router (no data files in context, inspect-before-ask, reference sanity gate, `pipefail` + per-step validation, pin and record versions, estimate before spending, thresholds are defaults, controlled-access data).
- Shared scripts: `peek.sh` (bounded preview of FASTQ/BAM/CRAM/VCF/BED/GTF), `ref_sanity.sh` (contig name/length/build consistency gate), `check_env.sh` (tool/version/scheduler probe used via dynamic context injection).
- `paths:` frontmatter so skills activate on relevant files; `genomics` owns the generic globs.
- CI: `scripts/lint_skills.py` (frontmatter, size, bash-block syntax via shellcheck, deprecated tools, nf-core registry, internal paths), `scripts/test_scripts.sh`, lychee link check. Removed `settings.json`; added `LICENSE` (MIT, as the README already claimed).

### New skills
`genomics-longread`, `genomics-somatic`, `genomics-methylation`, `genomics-metagenomics`, `genomics-benchmark`, `genomics-workflow`.

### Fixed (commands that were broken or wrong)
- VEP command broken by comments inside backslash continuation; non-existent `--plugin ClinVar`; `--plugin gnomADc` misdescribed as constraint; filtering on `INFO/SYMBOL`/`INFO/gnomAD_AF` (fields live in CSQ → `bcftools +split-vep`); `INFO/FILTER`.
- ATAC-seq: MACS `-f BAMPE` combined with `--shift/--extsize` on already Tn5-shifted reads; gawk-only `and()` read shift that broke mate coordinates; `grep -v chrM` on SAM text; standalone HMMRATAC; `TSSEscore()` given a path.
- GATK gCNV/somatic CNV missing `PreprocessIntervals`; `FilterIntervals --minimum-mappability` without a mappability track.
- Mutect2 shown without orientation-bias model, contamination estimation or PoN construction.
- DeepVariant run on BQSR output; `featureCounts` without `-s`; Salmon index without decoys; `lfcShrink(coef=)` without `relevel()`; DiffBind `minMembers` trap.
- `mosdepth -b … --by 1`; `multiqc -m featureCounts`; `verifybamid2 --Ref`; `nf-core/fastqc` and `nf-core/copynumber` (do not exist); R `copynumber` (removed from Bioconductor); `wget` HTTPS glob; NCBI dbSNP contig names; sarek `--genome` inconsistency and `--wes` without `--intervals`.

### Updated
DeepVariant 1.10, VEP/Ensembl 116, GENCODE v50 / vM39, gnomAD v4.1 (730,947 exomes + 76,215 genomes), MACS3, GATK 4.7 with NVScoreVariants, JASPAR2024, GRCm39, ClinGen SVI calibrated PP3/BP4 thresholds (REVEL, CADD, SpliceAI, AlphaMissense), ACMG without PP5/BP6, smoove for LUMPY, slivar for GEMINI, ENCODE TSS-enrichment bands.

## 1.0.1 — 2026-01-22
- Fix reference file URLs to use verified public sources.

## 1.0.0 — 2026-01-22
- Initial release: genomics, wgs, rnaseq, chipseq, annotation, qc, cnv.
