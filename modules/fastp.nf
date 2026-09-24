process FASTP {
    conda "${params.unified_conda_env}"
    tag "fastp on $sample_id"
    label = 'io_net'
    publishDir "${params.outdir}/1_fastp/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'

    publishDir "${params.outdir}/1_fastp/log_files/${sample_id}", 
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.clean.1.fastq.gz"), path("${sample_id}.clean.2.fastq.gz"), emit: reads
    path '.command.*'

    script:
    """
    # module load fastp/0.23.4-GCC-13.2.0
    fastp -w 4 -q 25 \
	--cut_front --cut_tail \
	-l 20 --detect_adapter_for_pe \
	-i $fq1 -I $fq2 \
	-o ${sample_id}.clean.1.fastq.gz \
	-O ${sample_id}.clean.2.fastq.gz \
	-h ${sample_id}_fastp.html \
	-j ${sample_id}_fastp.json \
	--thread 4
    """
}
