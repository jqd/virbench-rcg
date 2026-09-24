process GATHER_RESULTS {
    tag "$sample_id - $tool"
    label = 'io_limited'

    publishDir "${params.outdir}/13_combined_results/${tool}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), 
          path(blast_file),
          val(tool),
          path(kaiju_file),
          path(mmseqs2_file),
          path(kraken2_file), 
          path(fve_file), 
          path(viromescan_file)

    output:
    tuple val(sample_id), path("${sample_id}_*"), val(tool), emit: results

    script:
    """
    # Create copies with sample_id prefix for clarity
    cp ${blast_file} ${sample_id}_blast.tsv
    cp ${kaiju_file} ${sample_id}_kaiju.tsv
    cp ${mmseqs2_file} ${sample_id}_mmseqs2.tsv
    cp ${kraken2_file} ${sample_id}_kraken2.tsv
    cp ${fve_file} ${sample_id}_fve.tsv
    cp ${viromescan_file} ${sample_id}_viromescan.tsv
    """
}
