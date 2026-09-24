process MAP_BACK_TO_CONTIGS {
    conda "${params.unified_conda_env}"
    tag "bowtie2 on $sample_id"
    label = 'io_mem'

    publishDir "${params.outdir}/5_contig_abun/${tool}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.tsv'

    publishDir "${params.outdir}/5_contig_abun/log_files/${tool}/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}



    input:
    tuple val(sample_id), path(fq1), path(fq2), path(contigs), val(tool)

    output:
    tuple val(sample_id), path(contigs), path("${sample_id}_${tool}.abundance.tsv"), val(tool), emit: abundance
    path '.command.*'

    script:
    """
    # module load bowtie2/2.5.3
    # module load samtools/1.17
    bowtie2-build $contigs contigs_index
    bowtie2 -x contigs_index -1 $fq1 -2 $fq2 -S ${sample_id}.sam --very-sensitive -p 4
    samtools view -bS ${sample_id}.sam > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}_sorted_pos.bam
    samtools index ${sample_id}_sorted_pos.bam
    samtools idxstats ${sample_id}_sorted_pos.bam | cut -f1,3 | grep -v "^*" > ${sample_id}_${tool}.abundance.tsv
    """
}

