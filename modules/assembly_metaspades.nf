process ASSEMBLE_METASPADES {
    conda "${params.unified_conda_env}"
    tag "metaSpades on $sample_id"
    label = 'io_mem'
    publishDir "${params.outdir}/4_metaspades/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.scaffolds.fasta'

    publishDir "${params.outdir}/4_metaspades/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.scaffolds.fasta"), emit: contigs
    path '.command.*'

    script:
   // tool = getSoftwareName(task.process)
    """
    metaspades.py -1 $fq1 -2 $fq2 -o . -t 16
    cp scaffolds.fasta ${sample_id}.scaffolds.fasta
    """
}

