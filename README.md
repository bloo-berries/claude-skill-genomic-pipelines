# Genomics skills for Claude Code

Thirteen [Claude Code skills](https://code.claude.com/docs/en/skills) that turn Claude into a careful bioinformatics engineer for sequencing data: short- and long-read variant calling, somatic analysis, RNA-seq (bulk, single-cell, spatial), ChIP/ATAC/CUT&RUN, CNV/SV, annotation and interpretation, methylation, metagenomics, QC, benchmarking and workflow engineering.

The skills are written for an **autonomous agent**, not for a human reading a wiki: every skill starts from a fixed set of operating rules (never load data files into context, inspect headers before asking, verify reference consistency before combining inputs, validate each step, pin and record versions, treat thresholds as defaults, treat human data as controlled-access) and ships bounded inspection scripts so Claude can look at a 200 GB directory without reading it.

## Install

Project-wide (committed with the repository, shared by the team):

```bash
git clone https://github.com/bloo-berries/claude-skill-genomic-pipelines.git
cp -r claude-skill-genomic-pipelines/.claude/skills /path/to/your/project/.claude/
```

Personal (every project on this machine):

```bash
cp -r claude-skill-genomic-pipelines/.claude/skills/* ~/.claude/skills/
```

Restart Claude Code (or wait for live reload) and type `/` to see the skills. No `settings.json` entry is required; Claude Code discovers `.claude/skills/<name>/SKILL.md` automatically.

## Skills

| Skill | Invoke | Scope |
|---|---|---|
| Router | `/genomics` | Operating rules, inspection scripts, routing, shared conventions, reference/annotation/database sources (`reference.md`) |
| WGS/WES germline | `/genomics-wgs` | BWA-MEM2, GATK 4.7 (HaplotypeCaller, VQSR, NVScoreVariants, hard filters), DeepVariant 1.10, GLnexus, DRAGEN-GATK, sarek |
| Somatic | `/genomics-somatic` | Mutect2 full best practice, Strelka2, consensus, PoN, ASCAT/FACETS/PURPLE, MSI/TMB/HRD, signatures, AMP/ASCO/CAP tiers, PCGR, UMI consensus |
| RNA-seq | `/genomics-rnaseq` | Strandedness inference, decoy-aware Salmon, STAR, DESeq2/edgeR/limma, batch handling, splicing, fusions; scRNA (Seurat v5, Scanpy, ambient/doublet, integration, pseudobulk DE), spatial |
| ChIP / ATAC / CUT&RUN | `/genomics-chipseq` | MACS3, IDR, deepTools, DiffBind, motifs; ATAC Tn5 handling done correctly, `macs3 hmmratac`, TOBIAS, scATAC; SEACR and spike-in scaling |
| Annotation & interpretation | `/genomics-annotation` | VEP 116 with plugins, `bcftools +split-vep`, SnpEff, ANNOVAR, liftover; ACMG/ClinGen SVI with calibrated thresholds; slivar, Exomiser, STR, PGx, nf-core/raredisease |
| QC | `/genomics-qc` | fastp/Trim Galore/cutadapt, FastQC/MultiQC, mosdepth, contamination (VerifyBamID2, FastQ Screen, Kraken2), identity/sex/relatedness (somalier, peddy), PASS/WARN/FAIL tables |
| CNV / SV | `/genomics-cnv` | GATK gCNV and somatic CNV (with `PreprocessIntervals`), CNVkit, ExomeDepth/DECoN, Manta, smoove, Tiddit, dysgu, merging, AnnotSV/SVAFotate, ClinGen CNV scoring |
| Long read | `/genomics-longread` | dorado, pbmm2/minimap2, Clair3/DeepVariant, WhatsHap/HiPhase, Sniffles2/cuteSV/pbsv, modkit, TRGT/Straglr, hifiasm/Flye, Iso-Seq/bambu |
| Methylation | `/genomics-methylation` | Bismark, bwa-meth + MethylDackel, conversion QC, DSS/methylKit/dmrseq, arrays, nf-core/methylseq |
| Metagenomics | `/genomics-metagenomics` | Host removal, Kraken2/Bracken, MetaPhlAn 4, MAGs (MEGAHIT → binning → CheckM2/GTDB-Tk), HUMAnN, AMR, 16S/QIIME 2, nf-core/taxprofiler, mag, ampliseq |
| Benchmarking | `/genomics-benchmark` | GIAB truth sets and stratifications, hap.py / vcfeval, Truvari, SEQC2 somatic, nf-core/variantbenchmarking, minimum F1 by assay |
| Workflow engineering | `/genomics-workflow` | Nextflow DSL2 and nf-core mechanics, Snakemake 8, SLURM, containers and env pinning, cloud and accelerated stacks, failure triage, templates |

Skills also activate automatically when Claude works on matching files (`*.fastq.gz`, `*.bam`, `*.vcf.gz`, `*.narrowPeak`, `*.h5ad`, `*.pod5`, `*.nf`, `Snakefile`, …) or when the request matches a skill's description.

## Layout

```text
.claude/skills/
├── genomics/            SKILL.md  reference.md  scripts/{peek.sh,ref_sanity.sh,check_env.sh}
├── genomics-wgs/        SKILL.md  reference.md
├── genomics-somatic/    SKILL.md  reference.md
├── genomics-rnaseq/     SKILL.md  bulk.md  singlecell.md
├── genomics-chipseq/    SKILL.md  atac.md  cutandrun.md
├── genomics-annotation/ SKILL.md  acmg.md  rare-disease.md
├── genomics-qc/         SKILL.md
├── genomics-cnv/        SKILL.md  reference.md
├── genomics-longread/   SKILL.md  reference.md
├── genomics-methylation/SKILL.md  dmr.md
├── genomics-metagenomics/SKILL.md
├── genomics-benchmark/  SKILL.md
└── genomics-workflow/   SKILL.md  templates.md
scripts/lint_skills.py   repo lint (run in CI)
scripts/test_scripts.sh  smoke test for the shared scripts
```

Each `SKILL.md` has the same sections: operating rules, environment probe (live tool/version detection), inspect-first table, decision table, workflow with per-step validation, QC thresholds with rationale, pitfalls, nf-core pipelines (names validated in CI), validation, additional resources.

## Research use only

Variant classifications, CNV scoring and somatic tiers produced with these skills are research outputs. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189 or local equivalent) with the exact ClinVar/gnomAD/transcript-set versions recorded. Treat all human sequencing data as controlled-access; the skills instruct Claude not to upload data or paste identifiers, but the data-use agreement is yours to enforce.

## Verified against (2026-11-22)

GATK 4.7.0.0, DeepVariant 1.10.0, Ensembl/VEP 116, GENCODE v50 / vM39, gnomAD v4.1, MACS3, samtools/bcftools 1.2x, nf-core/sarek 3.12, Nextflow 25.x, Snakemake 8/9. Versions drift; the monthly CI run validates nf-core names and links, and every skill tells Claude to run `tool --version` rather than trust a number in a file.

## Development

```bash
python3 scripts/lint_skills.py      # frontmatter, size, bash syntax (needs shellcheck), deprecated tools, nf-core registry, paths
bash scripts/test_scripts.sh        # shared-script smoke test
```

Pull requests: keep `SKILL.md` under 400 lines, move long command dumps into supporting files, verify every version string against a primary source, and never recommend a tool the lint lists as deprecated. See `CHANGELOG.md`.

## License

MIT — see `LICENSE`.
