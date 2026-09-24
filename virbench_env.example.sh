#!/usr/bin/env bash
# Copy to virbench_env.sh, edit paths, then run:
#   source ./virbench_env.sh

# Shared references
export VIRBENCH_UNIFIED_CONDA_ENV="$HOME/.conda/envs/virbench-rcg"
export VIRBENCH_FVE_DIR="$HOME/bin/FastViromeExplorer"
export VIRBENCH_VIROMESCAN_CONDA_ENV="$HOME/.conda/envs/viromescan"
export VIRBENCH_HOST_INDEX="/path/to/host/bowtie2/index/prefix"
export VIRBENCH_NCBI_TAXONOMY_DIR="/path/to/taxonomy"
export VIRBENCH_GROUND_TRUTH="/path/to/ground_truth.tsv"

# Default DB profile paths
export VIRBENCH_DEFAULT_DB_ROOT="/path/to/default_db"
export VIRBENCH_DEFAULT_BLAST_DB="${VIRBENCH_DEFAULT_DB_ROOT}/nt/nt"
export VIRBENCH_DEFAULT_KAIJU_DB="${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25/kaiju_db_nr.fmi"
export VIRBENCH_DEFAULT_KAIJU_TAX_PATH="${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25"
export VIRBENCH_DEFAULT_MMSEQS2_DB="${VIRBENCH_DEFAULT_DB_ROOT}/mmseqs_NR/mmseqs_NR"
export VIRBENCH_DEFAULT_FVE_DB="${VIRBENCH_DEFAULT_DB_ROOT}/fve/ncbi-virus-kallisto-index-k31.idx"
export VIRBENCH_DEFAULT_FVE_VLIST="${VIRBENCH_DEFAULT_DB_ROOT}/fve/ncbi-viruses-list.txt"
export VIRBENCH_DEFAULT_FVE_ACC2TAXON="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data/fastviromeexplorer/ncbi-viruses-acc_8957_to_taxon.tsv"
export VIRBENCH_DEFAULT_VIROMESCAN_DB="virus_ALL"
export VIRBENCH_DEFAULT_VIROMESCAN_ACC2TAXON="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data/viromescan/VirusALLcomplete_4370_acc_to_taxon.tsv"
export VIRBENCH_DEFAULT_KRAKEN2_DB="${VIRBENCH_DEFAULT_DB_ROOT}/kraken2_db/k2_standard_20251015"

# RefSeq-CG profile paths
export VIRBENCH_REFSEQ_DB_ROOT="/path/to/60774_RefSeq_CG_db"
export VIRBENCH_REFSEQ_BLAST_DB="${VIRBENCH_REFSEQ_DB_ROOT}/blastn_refseq_CG/blast_60774_RefSeq_CG"
export VIRBENCH_REFSEQ_KAIJU_DB="${VIRBENCH_REFSEQ_DB_ROOT}/kaiju_refseq_CG/kaiju_db.refseq_60774.bwt.fmi"
export VIRBENCH_REFSEQ_KAIJU_TAX_PATH="${VIRBENCH_REFSEQ_DB_ROOT}/kaiju_refseq_CG/taxonomy"
export VIRBENCH_REFSEQ_MMSEQS2_DB="${VIRBENCH_REFSEQ_DB_ROOT}/mmseqs2_refseq_CG/60774_refseq_mmseqs_db"
export VIRBENCH_REFSEQ_FVE_DB="${VIRBENCH_REFSEQ_DB_ROOT}/FVE_refseq_CG/60774_RefSeq_CG_trans_kallisto_index.idx"
export VIRBENCH_REFSEQ_FVE_VLIST="${VIRBENCH_REFSEQ_DB_ROOT}/FVE_refseq_CG/60774_RefSeq_CG_trans_list_renamed.tsv"
export VIRBENCH_REFSEQ_VIROMESCAN_DB="60774_RefSeq_CG"
export VIRBENCH_REFSEQ_KRAKEN2_DB="${VIRBENCH_REFSEQ_DB_ROOT}/kraken2_refseq_CG"
