process ASSEMBLE_MEGAHIT {
    conda "${params.unified_conda_env}"
    tag "megahit on $sample_id"
    label = 'mem_medium'
    publishDir "${params.outdir}/3_megahit/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.contigs.fa'

    publishDir "${params.outdir}/3_megahit/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}
    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.fa"), emit: contigs 
    path '.command.*'

    script:
    """
    megahit -t 4 --force \
        -1 ${fq1} \
        -2 ${fq2} \
        --out-prefix ${sample_id} \
        -o .
    """
}

