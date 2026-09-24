//k2 dones't support gzipped input files
//use kraken2 to annotate reads


process KRAKEN2_ANNOTATE {
    conda "${params.unified_conda_env}"
    tag "Kraken2 on $sample_id"
    label 'io_mem'

    publishDir "${params.outdir}/11_kraken2/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}.${tool}*"

    publishDir "${params.outdir}/11_kraken2/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2) 

    output:
    tuple val(sample_id), path("${sample_id}.${tool}_abundance.tsv"), emit: results
    path("${sample_id}.${tool}*")
    path '.command.*'


    script:
    tool = "kraken2"
    """
    # module load kraken/2.1.2
    export KRAKEN_NUM_THREADS=${task.cpus}
    unset OMP_NUM_THREADS
    export OMP_NUM_THREADS=${task.cpus}
     
    kraken2 --db ${params.kraken2_db} \
          --paired --threads ${task.cpus} --report ${sample_id}.${tool}_taxonomy.txt \
          --output ${sample_id}.${tool}_output.txt \
          $fq1 $fq2   

    # Generate abundance table from kraken2 output
    (echo -e "Taxon\\tAbundance"; \
     cut -f3 ${sample_id}.${tool}_output.txt | \
        sort | uniq -c | \
        awk 'BEGIN{OFS="\\t"} {print \$2, \$1*2}') \
        > ${sample_id}.${tool}_abundance.tsv
    """
}

// C       M05128:454:000000000-D8875:1:1101:27814:11591   645687  151|151 645687:12 39733:5 645687:100 |:| 645687:117
// C       M05128:454:000000000-D8875:1:1102:25819:13486   645687  136|136 645687:102 |:| 645687:102
// C       M05128:454:000000000-D8875:1:1102:15605:11445   645687  151|151 645687:117 |:| 645687:57 0:39 645687:21
// C       M05128:454:000000000-D8875:1:1102:18230:19645   645687  137|137 645687:52 0:39 645687:12 |:| 645687:12 0:39 645687:52