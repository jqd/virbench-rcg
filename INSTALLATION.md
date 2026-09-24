# VirBench-RCG installation notes

## 1. Create the shared environment

From the repository root:

```bash
conda env create -f envs/virbench-rcg.yml
```

## 2. Install FastViromeExplorer

Follow the upstream installation instructions:

<https://github.com/saima-tithi/FastViromeExplorer>

Install FastViromeExplorer locally, for example:

```bash
cd "$HOME/bin"
git clone https://github.com/saima-tithi/FastViromeExplorer.git FastViromeExplorer
cd FastViromeExplorer
```

Add the FastViromeExplorer installation directory to `virbench_env.sh`:

```bash
export VIRBENCH_FVE_DIR="$HOME/bin/FastViromeExplorer"
```

The pipeline uses this directory as `${FVE_dir}/bin` for the Java classpath.
The same setting is included in `virbench_env.example.sh`; source the file
after editing your paths.

### Kallisto 0.43.1

FastViromeExplorer requires Kallisto 0.43.1. The FastViromeExplorer repository provides the Linux executable in `tools-linux`. Copy it into the shared `virbench-rcg` environment:

```bash
cp "$HOME/bin/FastViromeExplorer/tools-linux/kallisto" \
   "$HOME/.conda/envs/virbench-rcg/bin/"
chmod +x "$HOME/.conda/envs/virbench-rcg/bin/kallisto"
"$HOME/.conda/envs/virbench-rcg/bin/kallisto"
```

The final command should report version `0.43.1`.

## 3. Install ViromeScan

Follow the upstream installation instructions:

<https://github.com/simonerampelli/viromescan>

The pipeline uses the environment at:

```text
~/.conda/envs/viromescan
```

For a custom ViromeScan database, place the Bowtie2 database under:

```text
~/.conda/envs/viromescan/viromescan/database/bowtie2/
```

The default ViromeScan accession-to-taxonomy mapping file is included in this
repository at `data/viromescan/VirusALLcomplete_4370_acc_to_taxon.tsv`; users do
not need to download it separately. `virbench_env.example.sh` configures this
path automatically.

After installation, edit:

```text
~/.conda/envs/viromescan/viromescan/viromescan.sh
```

At the end of that file, comment out these post-processing lines so ViromeScan does not overwrite or remove pipeline output directories:

```bash
# Rscript ../viromescan.R
# cd ../../
# mv $OUTPUT_DIR/results .
# rm -rf $OUTPUT_DIR/*
# mv results/* $OUTPUT_DIR/
# rm -rf results
```

## 4. Set up default databases (db_default)

This section describes how to prepare the databases used by the `db_default`
profile and how to map them to `VIRBENCH_DEFAULT_*` variables in
`virbench_env.sh`.

Working example:

```bash
cp virbench_env.example.sh virbench_env.sh
```

Edit `virbench_env.sh` and set all `VIRBENCH_DEFAULT_*` paths.

Use a base variable:

```bash
export VIRBENCH_DEFAULT_DB_ROOT="/path/to/default_db"
export VIRBENCH_DEFAULT_BLAST_DB="${VIRBENCH_DEFAULT_DB_ROOT}/nt/<db_prefix>"
export VIRBENCH_DEFAULT_KAIJU_DB="${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25/kaiju_db_nr.fmi"
export VIRBENCH_DEFAULT_KAIJU_TAX_PATH="${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25"
export VIRBENCH_DEFAULT_MMSEQS2_DB="${VIRBENCH_DEFAULT_DB_ROOT}/mmseqs_NR/<mmseqs_db_prefix>"
export VIRBENCH_DEFAULT_FVE_DB="${VIRBENCH_DEFAULT_DB_ROOT}/fve/ncbi-virus-kallisto-index-k31.idx"
export VIRBENCH_DEFAULT_FVE_VLIST="${VIRBENCH_DEFAULT_DB_ROOT}/fve/ncbi-viruses-list.txt"
export VIRBENCH_DEFAULT_VIROMESCAN_DB="virus_ALL"
export VIRBENCH_DEFAULT_KRAKEN2_DB="${VIRBENCH_DEFAULT_DB_ROOT}/kraken2_db/<kraken2_db_dir>"
```

The accession-to-taxonomy mappings for FVE and ViromeScan are configured
separately in `virbench_env.example.sh` and can remain repo-relative.

### 4.1 BLAST (blastn 2.12.0)

- Study DB: `Default BLAST DB`
- Link: <https://ftp.ncbi.nlm.nih.gov/blast/db/>
- Env var: `VIRBENCH_DEFAULT_BLAST_DB`

Download and prepare a BLAST nucleotide DB prefix:

```bash
mkdir -p "${VIRBENCH_DEFAULT_DB_ROOT}/nt"
cd "${VIRBENCH_DEFAULT_DB_ROOT}/nt"

# Example: download preformatted BLAST DB archives from NCBI FTP.
# Choose files matching your target DB and unpack them.

# Set prefix path in virbench_env.sh (example)
# export VIRBENCH_DEFAULT_BLAST_DB="${VIRBENCH_DEFAULT_DB_ROOT}/nt/<db_prefix>"
```

`VIRBENCH_DEFAULT_BLAST_DB` must point to the BLAST database prefix used by
`blastn` (not a directory alone).

### 4.2 Kaiju (1.7.3)

- Study DB: `Default Kaiju DB`
- Link: <https://bioinformatics-centre.github.io/kaiju/downloads.html>
- Env vars: `VIRBENCH_DEFAULT_KAIJU_DB`, `VIRBENCH_DEFAULT_KAIJU_TAX_PATH`

Download Kaiju DB and taxonomy:

```bash
mkdir -p "${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25"
cd "${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25"

# Example archive from the study table:
# kaiju_db_nr_2024-08-25.tgz
# tar -xzf kaiju_db_nr_2024-08-25.tgz

# Set in virbench_env.sh (examples)
# export VIRBENCH_DEFAULT_KAIJU_DB="${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25/kaiju_db_nr.fmi"
# export VIRBENCH_DEFAULT_KAIJU_TAX_PATH="${VIRBENCH_DEFAULT_DB_ROOT}/kaiju_db_nr_2024-08-25"
```

`VIRBENCH_DEFAULT_KAIJU_DB` should be the `.fmi` file path.
`VIRBENCH_DEFAULT_KAIJU_TAX_PATH` should be the directory containing taxonomy
files required by Kaiju.

### 4.3 MMseqs2 taxonomy (15.6f452)

- Study DB: `Default MMseqs2 DB`
- Link: <https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA>
- Env var: `VIRBENCH_DEFAULT_MMSEQS2_DB`

Prepare protein FASTA inputs and build MMseqs2 taxonomy DB:

```bash
mkdir -p "${VIRBENCH_DEFAULT_DB_ROOT}/mmseqs_NR"
cd "${VIRBENCH_DEFAULT_DB_ROOT}/mmseqs_NR"

# Download protein FASTA sources from NCBI FTP, then build DB with mmseqs2.
# Set in virbench_env.sh (example)
# export VIRBENCH_DEFAULT_MMSEQS2_DB="${VIRBENCH_DEFAULT_DB_ROOT}/mmseqs_NR/<mmseqs_db_prefix>"
```

`VIRBENCH_DEFAULT_MMSEQS2_DB` must point to the MMseqs2 database prefix used in
the pipeline.

### 4.4 Kraken2 (2.1.2)

- Study DB: `Default Kraken2 DB`
- Link: <https://benlangmead.github.io/aws-indexes/k2>
- Env var: `VIRBENCH_DEFAULT_KRAKEN2_DB`

Download and unpack the Kraken2 index:

```bash
mkdir -p "${VIRBENCH_DEFAULT_DB_ROOT}/kraken2"
cd "${VIRBENCH_DEFAULT_DB_ROOT}/kraken2"

# Example file from the study table:
# k2_standard_20251015.tar.gz
# tar -xzf k2_standard_20251015.tar.gz

# Set in virbench_env.sh (example)
# export VIRBENCH_DEFAULT_KRAKEN2_DB="${VIRBENCH_DEFAULT_DB_ROOT}/kraken2/<kraken2_db_dir>"
```

`VIRBENCH_DEFAULT_KRAKEN2_DB` should be the DB directory consumed by `kraken2`.

### 4.5 FastViromeExplorer (1.3)

- Study DB: `Default FVE index`
- Link: <https://code.vt.edu/saima5/FastViromeExplorer>
- Env vars: `VIRBENCH_DEFAULT_FVE_DB`, `VIRBENCH_DEFAULT_FVE_VLIST`

Prepare the Kallisto index and virus list:

```bash
mkdir -p "${VIRBENCH_DEFAULT_DB_ROOT}/fve"
cd "${VIRBENCH_DEFAULT_DB_ROOT}/fve"

# Example file from the study table:
# ncbi-virus-kallisto-index-k31.idx

# Set in virbench_env.sh (examples)
# export VIRBENCH_DEFAULT_FVE_DB="${VIRBENCH_DEFAULT_DB_ROOT}/fve/ncbi-virus-kallisto-index-k31.idx"
# export VIRBENCH_DEFAULT_FVE_VLIST="${VIRBENCH_DEFAULT_DB_ROOT}/fve/<virus_list.tsv>"
```

### 4.6 ViromeScan

- Study DB: `Default ViromeScan DB`
- Link: <https://sourceforge.net/projects/viromescan/>
- Env var: `VIRBENCH_DEFAULT_VIROMESCAN_DB`

For default-db runs in this pipeline, set:

```bash
export VIRBENCH_DEFAULT_VIROMESCAN_DB="virus_ALL"
```

This value is already provided in `virbench_env.example.sh`.


## 5. Set up RefSeq-CG databases (db_refseq)

RefSeq-CG databases for this pipeline can be downloaded from Harvard Dataverse:

<https://doi.org/10.7910/DVN/BQG680>

Because of Dataverse size limits, each database is split into multiple
small `.tar.gz` pieces. For each database, download all pieces first, then
merge them into one archive and extract it.

### 5.1 Merge split archives and extract

For each database folder downloaded from Dataverse:

```bash
# Example working directory containing split files for one database
cd /path/to/downloaded_parts/<database_name>

# 1) Verify all pieces are present
ls -1

# 2) Merge split pieces in version order (adjust pattern if needed)
printf '%s\n' <database_name>.tar.gz.part* | sort -V | xargs cat -- > <database_name>.tar.gz

# 3) Extract
tar -xzf <database_name>.tar.gz
```

Repeat this process for BLAST, Kaiju, MMseqs2, FVE, Kraken2, and ViromeScan
until all database bundles have been reconstructed and extracted.

After extraction is complete, choose the parent directory that contains the
RefSeq-CG database tree, then set `VIRBENCH_REFSEQ_*` variables in
`virbench_env.sh` using the structure below.

### 5.2 Example layout and exports

Assume the extracted directory is:

```text
/path/to/60774_RefSeq_CG_db
```

Set variables using the base path:

```bash
export VIRBENCH_REFSEQ_DB_ROOT="/path/to/60774_RefSeq_CG_db"
export VIRBENCH_REFSEQ_BLAST_DB="${VIRBENCH_REFSEQ_DB_ROOT}/blastn_refseq_CG/blast_60774_RefSeq_CG"
export VIRBENCH_REFSEQ_KAIJU_DB="${VIRBENCH_REFSEQ_DB_ROOT}/kaiju_refseq_CG/kaiju_db.refseq_60774.bwt.fmi"
export VIRBENCH_REFSEQ_KAIJU_TAX_PATH="${VIRBENCH_REFSEQ_DB_ROOT}/kaiju_refseq_CG/taxonomy"
export VIRBENCH_REFSEQ_MMSEQS2_DB="${VIRBENCH_REFSEQ_DB_ROOT}/mmseqs2_refseq_CG/60774_refseq_mmseqs_db"
export VIRBENCH_REFSEQ_FVE_DB="${VIRBENCH_REFSEQ_DB_ROOT}/FVE_refseq_CG/60774_RefSeq_CG_trans_kallisto_index.idx"
export VIRBENCH_REFSEQ_FVE_VLIST="${VIRBENCH_REFSEQ_DB_ROOT}/FVE_refseq_CG/60774_RefSeq_CG_trans_list_renamed.tsv"
export VIRBENCH_REFSEQ_VIROMESCAN_DB="60774_RefSeq_CG"
export VIRBENCH_REFSEQ_KRAKEN2_DB="${VIRBENCH_REFSEQ_DB_ROOT}/kraken2_refseq_CG"
```

