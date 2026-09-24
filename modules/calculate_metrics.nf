process CALCULATE_METRICS {
    conda "${params.unified_conda_env}"
    tag "$sample_id - $tool"
    label = 'io_limited'

    publishDir "${params.outdir}/12_metrics/${tool}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.{pdf,tsv}'

    publishDir "${params.outdir}/12_metrics/${tool}/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), 
          path(blast_results),
          val(tool),
          path(kaiju_results), 
          path(mmseqs2_results), 
          path(kraken2_results), 
          path(fve_results), 
          path(viromescan_results),
          path(ground_truth)

    output:
    tuple val(sample_id), path("${sample_id}_metrics.tsv"), val(tool), emit: metrics
    tuple val(sample_id), path("${sample_id}_metrics_plot.pdf"), val(tool), emit: metrics_plot
    tuple val(sample_id), path("species_results_*.tsv"), val(tool), emit: table
    path '.command.*'


    script:
    """
    # module load R/4.4.0
    
    Rscript ${projectDir}/bin/calculate_metrics.R \
        --blast ${blast_results} \
        --kaiju ${kaiju_results} \
        --mmseqs2 ${mmseqs2_results} \
        --kraken2 ${kraken2_results} \
        --fve ${fve_results} \
        --viromescan ${viromescan_results} \
        --ground_truth ${ground_truth} \
        --ncbi_taxonomy_dir ${params.ncbi_taxonomy_dir} \
        --tool ${tool} \
        --output ${sample_id}_metrics.tsv \
        --FVE_default_db ${params.FVE_default_db} \
        --FVE_acc2taxon "${params.FVE_acc2taxon}" \
        --viromescan_default_db ${params.viromescan_default_db} \
        --viromescan_acc2taxon "${params.viromescan_acc2taxon}"
    """
}
