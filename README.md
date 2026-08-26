# Snakemake Reimplementation: CYP51 Resistance Mutation in Amphotericin B-Resistant *Leishmania mexicana*

A Snakemake pipeline that identifies genomic variants unique to an amphotericin B (AmpB)-resistant *Leishmania mexicana* strain, re-engineered from a manual NGS variant-calling workflow into a reproducible, automated pipeline.

## Background

This project re-runs the genomic variant-calling component of my earlier [Integrative Multi-Omics Analysis of Amphotericin B Resistance in *Leishmania mexicana*](https://github.com/Claire-bioinformatics/Integrative-Multi-Omics-Analysis-of-Amphotericin-B-Resistance-in-Leishmania-mexicana) project — originally performed as a series of manual command-line steps — as a fully automated Snakemake workflow. The original university HPC data is no longer accessible, so this version re-runs the pipeline end-to-end against public data, executed via Google Colab.

**Data source:** Mwenechanya et al. (2017), *Sterol 14α-demethylase mutation leads to amphotericin B resistance in Leishmania mexicana*, PLOS Neglected Tropical Diseases 11(6):e0005649. Raw WGS reads from ENA project [PRJEB10872](https://www.ebi.ac.uk/ena/browser/view/PRJEB10872) — runs `ERR1517168` (wild-type) and `ERR1517169` (AmpB-resistant), strain M379, Illumina GAIIx.

**Reference genome:** *L. mexicana* MHOM/GT/2001/U1103, NCBI assembly `GCA_000234665.4`.

## Pipeline overview

```
Raw reads (WT, AmpB)
   |
   |-- FastQC (quality control)
   |
   |-- Trim Galore (adapter trimming)
   |
   |-- Bowtie2 (alignment to reference genome)
   |
   |-- samtools (sort + index BAMs)
   |
   |-- bamaddrg + FreeBayes (joint variant calling)
   |
   |-- vcffilter (QUAL > 20, SNPs only)
   |
   |-- bcftools isec (identify variants unique to AmpB-resistant strain)
   |
   |-- SnpEff (variant annotation, custom L. mexicana database)
   |
   `-- SnpSift (extract non-synonymous / high-impact variants) -> results/snpeff/AmpB_unique_nonsyn.tsv
```

## Key finding

Of 52 non-synonymous/high-impact SNPs found uniquely in the AmpB-resistant strain, one lands directly in **CYP51** (sterol 14α-demethylase, gene `LMXM_11_1100`):

| Gene | Mutation | Effect | Genotype |
|---|---|---|---|
| CYP51 (LMXM_11_1100) | c.527A>T | p.Asn176Ile (N176I) | 1/1 (homozygous) |

This matches the resistance mechanism reported in the source publication — a mutation in the ergosterol biosynthesis pathway target enzyme, consistent with reduced AmpB binding affinity.

## Repository structure

```
├── Snakefile              # pipeline rules
├── config.yaml            # reference genome path, ploidy, SnpEff config path
├── config/
│   └── snpEff.config      # custom SnpEff database definition
├── snakemake.ipynb        # Colab notebook: environment setup + full pipeline run
└── README.md
```

## How to run

This pipeline was executed on Google Colab (free tier) due to its bioinformatics tool dependencies (conda-based) and moderate compute requirements. `snakemake.ipynb` handles the full setup and run:

1. Open `snakemake.ipynb` in Google Colab
2. Run cells sequentially — this installs a conda environment (Snakemake, FastQC, Trim Galore, Bowtie2, samtools, FreeBayes, vcflib, bcftools, SnpEff, SnpSift), clones this repo, downloads the public FASTQ + reference data, builds the custom SnpEff database, and runs the pipeline (`snakemake --cores 2`)
3. Results are written to `results/`, with the key output at `results/snpeff/AmpB_unique_nonsyn.tsv`

To run locally instead (with the tools installed via conda/mamba and data already downloaded into `data/`):

```bash
snakemake --cores <N>
```

## Tools

Snakemake · FastQC · Trim Galore · Bowtie2 · samtools · FreeBayes · vcflib · bcftools · SnpEff · SnpSift

## Future work

- **Biological replicates** — currently one WT vs. one AmpB-resistant sample; validating CYP51 N176I across multiple resistant isolates would confirm it's not strain-specific
- **Cross-validate variant calls** — run a second caller (e.g. GATK HaplotypeCaller or DeepVariant) alongside FreeBayes to check concordance on the key SNPs
- **Structural modeling of the resistance mutation** — predict the mutant CYP51 structure (AlphaFold2/ColabFold) and dock amphotericin B against WT vs. mutant to assess how N176I plausibly disrupts drug binding
- **Extend Snakemake beyond genomics** — Snakemake is a general-purpose workflow manager, not genomics-specific. The proteomics and lipidomics arms of the parent multi-omics project (currently manual analysis) could be brought into this same reproducible framework, orchestrating mass-spec search tools (e.g. MaxQuant) and structural prediction/docking steps alongside the existing variant-calling rules

