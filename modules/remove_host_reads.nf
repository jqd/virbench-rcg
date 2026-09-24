process REMOVE_HOST_READS {
    conda "${params.unified_conda_env}"
    tag "bowtie2 $sample_id"
    label = 'io_net'
    publishDir "${params.outdir}/2_bowtie2/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'

    publishDir "${params.outdir}/2_bowtie2/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}
    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}_host_removed_R1.fastq.gz"), path("${sample_id}_host_removed_R2.fastq.gz"), emit: reads
    path '.command.*'

    script:
    """
    bowtie2 -x ${params.host} \
     -1 $fq1 -2 $fq2 \
     --very-sensitive-local --threads 5 \
     --un-conc-gz ${sample_id}_host_removed \
     -S ${sample_id}_mapped_and_unmapped.sam

    mv ${sample_id}_host_removed.1 ${sample_id}_host_removed_R1.fastq.gz
    mv ${sample_id}_host_removed.2 ${sample_id}_host_removed_R2.fastq.gz

    """
}

