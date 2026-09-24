# VirBench-RCG

VirBench-RCG is a Nextflow DSL2 pipeline for viral metagenomics profiling and benchmarking. It combines read-based and contig-based taxonomic annotation and supports two database strategies:

- tool-specific default databases (`db_default`)
- a uniform RefSeq-CG database set (`db_refseq`)

## Table of contents

- [Workflow overview](#workflow-overview)
- [Requirements](#requirements)
- [Installation and setup](#installation-and-setup)
- [Profiles and execution modes](#profiles-and-execution-modes)
- [Input sample sheet](#input-sample-sheet)
- [Usage examples](#usage-examples)
- [Parameter reference (--help)](#parameter-reference---help)
- [Required vs optional parameters](#required-vs-optional-parameters)
- [Environment variable mapping](#environment-variable-mapping)
- [Output structure](#output-structure)
- [Validation checklist](#validation-checklist)
- [Repository structure](#repository-structure)
- [Citation](#citation)
- [Notes for GitHub release](#notes-for-github-release)

## Workflow overview

The pipeline performs:

- Read quality filtering with Fastp
- Optional host-read removal with Bowtie2
- Assembly with MEGAHIT, metaSPAdes, or both
- Read mapping back to contigs for abundance estimation
- CD-HIT contig clustering
- Contig annotation with BLAST, Kaiju, and MMseqs2 taxonomy
- Read-based profiling with FastViromeExplorer, ViromeScan, and Kraken2
- Taxonomy normalization with TaxonKit
- Optional benchmarking with precision, recall, F1 score, and L2 distance

## Requirements

- Nextflow 23.10.0 or newer
- Conda
- local installation of required databases and references

Create the shared Conda environment:

```bash
conda env create -f envs/virbench-rcg.yml
```

Detailed external-tool/database setup is in [INSTALLATION.md](INSTALLATION.md).

## Installation and setup

Copy and edit the environment file:

```bash
cp virbench_env.example.sh virbench_env.sh
# edit virbench_env.sh
source virbench_env.sh
```

The local `virbench_env.sh` is machine-specific and should not be committed.

Use [nextflow.config](nextflow.config) as the canonical configuration.

## Profiles and execution modes

Database profiles:

- `db_default`: tool-specific default databases
- `db_refseq`: uniform RefSeq-CG databases

Executor profiles (from [profiles/schedulers.config](profiles/schedulers.config)):

- `standard` or `local`: local execution
- `lsf`: LSF scheduler
- `slurm`, `pbs`, `sge`: available if your environment uses them

Typical profile combinations:

- local + default DB: `-profile standard,db_default`
- local + RefSeq DB: `-profile standard,db_refseq`
- LSF + default DB: `-profile lsf,db_default`
- LSF + RefSeq DB: `-profile lsf,db_refseq`

## Input sample sheet

Provide a tab-separated sample sheet with header:

```text
sample_id  fq1  fq2
```

`fq1` and `fq2` are paired-end FASTQ paths. See [samples.tsv](samples.tsv).

## Usage examples

Show built-in help:

```bash
nextflow -C nextflow.config run main.nf --help
```

Run with default databases:

```bash
nextflow -C nextflow.config run main.nf \
  -profile lsf,db_default \
  --samples samples.tsv \
  --assembler both \
  --skip_remove_host true \
  --outdir results_default
```

Run with RefSeq-CG databases:

```bash
nextflow -C nextflow.config run main.nf \
  -profile lsf,db_refseq \
  --samples samples.tsv \
  --assembler both \
  --outdir results_refseq
```

Run annotation-only (skip benchmarking metrics):

```bash
nextflow -C nextflow.config run main.nf \
  -profile lsf,db_default \
  --samples samples.tsv \
  --skip_metrics true \
  --outdir results_annotation_only
```

Resume interrupted runs with `-resume`.

## Parameter reference (--help)

The following list matches the `--help` message defined in [main.nf](main.nf).

Core options:

- `--samples <file>`: sample sheet TSV (default: `samples.tsv`)
- `--outdir <dir>`: output directory (help text default: `results`; profile defaults may override this)
- `--assembler <name>`: `megahit`, `metaspades`, or `both` (default: `both`)
- `--skip_fastp <true|false>`: skip read trimming (default: `false`)
- `--skip_metrics <true|false>`: skip metrics calculation (default: `false`)
- `--skip_remove_host <true|false>`: skip host removal (default: `false`)

Read preprocessing:

- `--host <bowtie2_index_prefix>`

Annotation database options:

- `--blast_db <path>`
- `--blast_evalue <value>` (default: `1e-5`)
- `--blast_number_bitscores <int>` (default: `5`)
- `--kaiju_db <path>`
- `--kaiju_tax_path <path>`
- `--mmseqs2_db <path>`
- `--kraken2_db <path>`
- `--FVE_db <path>`
- `--FVE_vlist <path>`
- `--FVE_dir <path>`
- `--FVE_default_db <true|false>`
- `--FVE_acc2taxon <path>`
- `--viromescan_db <name|path>`
- `--viromescan_default_db <true|false>`
- `--viromescan_acc2taxon <path>`

CD-HIT options:

- `--cdhit_identity <float>` (default: `0.95`)
- `--cdhit_word_length <int>` (default: `10`)

Taxonomy and metrics references:

- `--ncbi_taxonomy_dir <path>`
- `--taxonomy <path>` (optional precomputed lineage file)
- `--ground_truth <path>`

Other:

- `--help`

## Required vs optional parameters

Validation logic in [main.nf](main.nf) enforces the following:

Always required:

- `--outdir`
- `--host`
- `--blast_db`
- `--kaiju_db`
- `--kaiju_tax_path`
- `--mmseqs2_db`
- `--FVE_db`
- `--FVE_vlist`
- `--viromescan_db`
- `--kraken2_db`
- `--FVE_default_db`
- `--viromescan_default_db`

Required when metrics are enabled (`--skip_metrics false`):

- `--ncbi_taxonomy_dir`
- `--ground_truth`

Conditionally required:

- if `--FVE_default_db true`, then `--FVE_acc2taxon` is required
- if `--viromescan_default_db true`, then `--viromescan_acc2taxon` is required

## Environment variable mapping

Values in [virbench_env.sh](virbench_env.sh) populate Nextflow params through [nextflow.config](nextflow.config):

- `VIRBENCH_UNIFIED_CONDA_ENV` -> `params.unified_conda_env`
- `VIRBENCH_FVE_DIR` -> `params.FVE_dir`
- `VIRBENCH_VIROMESCAN_CONDA_ENV` -> `params.viromescan_conda_env`
- `VIRBENCH_HOST_INDEX` -> `params.host`
- `VIRBENCH_NCBI_TAXONOMY_DIR` -> `params.ncbi_taxonomy_dir`
- `VIRBENCH_GROUND_TRUTH` -> `params.ground_truth`

Default DB profile (`db_default`):

- `VIRBENCH_DEFAULT_BLAST_DB` -> `params.blast_db`
- `VIRBENCH_DEFAULT_KAIJU_DB` -> `params.kaiju_db`
- `VIRBENCH_DEFAULT_KAIJU_TAX_PATH` -> `params.kaiju_tax_path`
- `VIRBENCH_DEFAULT_MMSEQS2_DB` -> `params.mmseqs2_db`
- `VIRBENCH_DEFAULT_FVE_DB` -> `params.FVE_db`
- `VIRBENCH_DEFAULT_FVE_VLIST` -> `params.FVE_vlist`
- `VIRBENCH_DEFAULT_FVE_ACC2TAXON` -> `params.FVE_acc2taxon`
- `VIRBENCH_DEFAULT_VIROMESCAN_DB` -> `params.viromescan_db`
- `VIRBENCH_DEFAULT_VIROMESCAN_ACC2TAXON` -> `params.viromescan_acc2taxon`
- `VIRBENCH_DEFAULT_KRAKEN2_DB` -> `params.kraken2_db`

RefSeq profile (`db_refseq`):

- `VIRBENCH_REFSEQ_BLAST_DB` -> `params.blast_db`
- `VIRBENCH_REFSEQ_KAIJU_DB` -> `params.kaiju_db`
- `VIRBENCH_REFSEQ_KAIJU_TAX_PATH` -> `params.kaiju_tax_path`
- `VIRBENCH_REFSEQ_MMSEQS2_DB` -> `params.mmseqs2_db`
- `VIRBENCH_REFSEQ_FVE_DB` -> `params.FVE_db`
- `VIRBENCH_REFSEQ_FVE_VLIST` -> `params.FVE_vlist`
- `VIRBENCH_REFSEQ_VIROMESCAN_DB` -> `params.viromescan_db`
- `VIRBENCH_REFSEQ_KRAKEN2_DB` -> `params.kraken2_db`

## Output structure

Outputs are published under `--outdir` in tool/stage folders:

- `1_fastp`
- `2_bowtie2`
- `3_megahit`
- `4_metaspades`
- `5_contig_abun`
- `6_blast`
- `7_kaiju`
- `8_mmseqs2`
- `9_FVE`
- `10_viromescan`
- `11_kraken2`
- `12_metrics` (if metrics enabled)
- `13_combined_results`
- `cdhit_megahit` and/or `cdhit_metaspades`

Nextflow run reports are also generated in the project root:

- `report.html`
- `timeline.html`
- `trace.txt`

## Validation checklist

Before a full production run:

1. Source environment variables:
   `source virbench_env.sh`
2. Confirm help message:
   `nextflow -C nextflow.config run main.nf --help`
3. Run a one-sample smoke test:
  `nextflow -C nextflow.config run main.nf -profile lsf,db_default --samples samples.tsv --outdir smoke_test -resume`
4. Verify expected output folders and reports are created.

## Repository structure

- [main.nf](main.nf): workflow entry point
- [modules](modules): process modules
- [bin](bin): helper scripts
- [envs](envs): Conda environment specs
- [data](data): bundled small mapping/reference helper files
- [INSTALLATION.md](INSTALLATION.md): external tool/database setup

## Citation

If you use this workflow or associated benchmark resources, please cite:

Qidong Jia, Valerie Cortez, Amy Davis, Stacey Schultz-Cherry, Elisa B Margolis. 2026. "Comparative Assessment of Viral Profiling Tools for Known Virus Identification: Performance Gaps in Metagenomic Analysis." Harvard Dataverse. <https://doi.org/10.7910/DVN/BQG680>.

