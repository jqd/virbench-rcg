# Fetch Taxon ID from NCBI

A Python script to retrieve taxon IDs for NCBI accession numbers, with special handling for updated, obsolete, or replaced sequences.

## Features

- Fetches current taxon IDs for accession numbers using NCBI E-utilities
- Handles updated/replaced/obsolete sequences
- Detects when an accession has been replaced and reports the current accession
- Processes single accessions or batch files
- Respects NCBI API rate limits (max 3 requests/second without API key)
- Includes retry logic for robustness

## Requirements

```bash
# Python packages
pandas
```

Install with:
```bash
pip install pandas
```

Or use conda:
```bash
conda install pandas
```

## Usage

### Single Accession Query

```bash
python bin/fetch_taxid_from_ncbi.py \
    --accession NC_005817.1 \
    --email your.email@example.com
```

Output:
```
Query Accession: NC_005817.1
Current Accession: NC_005817.1
Taxon ID: 228391
Status: current
```

### Batch Processing from Simple List

Create a file `accessions.txt` with one accession per line:
```
NC_005817.1
AC_000007.1
AC_000008.1
NC_013060.1
```

Then run:
```bash
python bin/fetch_taxid_from_ncbi.py \
    --input accessions.txt \
    --email your.email@example.com \
    --output taxid_results.tsv
```

### Batch Processing from TSV File

If you have a TSV file with an accession column:
```bash
python bin/fetch_taxid_from_ncbi.py \
    --input blast_results.tsv \
    --column subject_accession \
    --email your.email@example.com \
    --output taxid_results.tsv
```

## Output Format

The script generates a TSV file with the following columns:

| Column | Description |
|--------|-------------|
| `query_accession` | The accession you queried |
| `taxid` | The taxon ID (if found) |
| `status` | Status: current, updated, replaced, not_found, invalid, or error |
| `current_accession` | The current accession (if different from query) |
| `error` | Error message (if applicable) |

### Status Codes

- **current**: Accession is current and active
- **updated**: Accession has been updated (version changed)
- **replaced**: Accession has been replaced by another
- **not_found**: Accession not found in NCBI database
- **invalid**: Invalid accession format
- **error**: Error occurred during retrieval

## Examples

### Example 1: Check if a sequence has been updated

```bash
python bin/fetch_taxid_from_ncbi.py \
    --accession NC_013060.1 \
    --email your@email.com
```

### Example 2: Process FVE results

Extract accessions from FVE output and fetch taxon IDs:
```bash
# Extract unique accessions from FVE results
cut -f1 FVE_sample1.abundance.tsv | tail -n +2 > fve_accessions.txt

# Fetch taxon IDs
python bin/fetch_taxid_from_ncbi.py \
    --input fve_accessions.txt \
    --email your@email.com \
    --output fve_taxids.tsv
```

### Example 3: Process BLAST results

```bash
# Extract unique subject accessions from BLAST output
grep -v '^#' sample1.blast.tsv | cut -f2 | sort -u > blast_accessions.txt

# Fetch taxon IDs
python bin/fetch_taxid_from_ncbi.py \
    --input blast_accessions.txt \
    --email your@email.com \
    --output blast_taxids.tsv \
    --delay 0.35
```

## Rate Limits

- **Without API key**: 3 requests per second
- **With API key**: 10 requests per second

Default delay is 0.4 seconds (2.5 requests/sec) to be safe. Adjust with `--delay` if needed.

To use an API key, modify the script to include:
```python
'api_key': 'YOUR_API_KEY'
```
in the URL parameters.

## Integration with Pipeline

You can integrate this into your Nextflow pipeline by creating a new process:

```groovy
process FETCH_NCBI_TAXIDS {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(accession_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_taxids.tsv")
    
    script:
    """
    python ${projectDir}/bin/fetch_taxid_from_ncbi.py \
        --input ${accession_file} \
        --email ${params.ncbi_email} \
        --output ${sample_id}_taxids.tsv
    """
}
```

## Troubleshooting

### "Too many requests" error
- Increase the `--delay` parameter
- Wait a few minutes before retrying

### "Accession not found"
- Check if the accession is correct
- The sequence might have been removed from NCBI
- Try searching on NCBI website manually

### Slow processing
- This is normal for large batches due to rate limiting
- Consider running overnight for thousands of accessions
- With NCBI API key, processing is ~3x faster

## Notes

- Always provide a valid email address (NCBI requirement)
- For large-scale processing, consider obtaining an NCBI API key
- The script includes automatic retry logic for transient failures
- Results are saved incrementally to prevent data loss on interruption
