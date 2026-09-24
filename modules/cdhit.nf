process CDHIT {
    conda "${params.unified_conda_env}"
    tag "CD-HIT on $sample_id"
    label = 'io_limited'

    publishDir path: "${params.outdir}", mode:'copy', overwrite: true, \
        pattern: "${outdir}/*"

    publishDir path: "${params.outdir}/${outdir}/log_files/", \
        mode:'copy', \
        overwrite: true, \
        pattern: '.command.*', \
        saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs), val(tool)

    output:
    tuple val(sample_id), path("${outdir}/${prefix}.clustered*.txt"), val(tool), emit: clusters
    tuple val(sample_id), path("${outdir}/${prefix}.clustered*.fa"), val(tool), emit: repseqs
    path '.command.*',                                    emit: command
   

    script:
        prefix = "cdhit_${sample_id}"
        outdir = "cdhit_${tool}/${prefix}"

    """
        mkdir -p $outdir
        cd-hit-est \
            -i $contigs \
            -c 0.8 -d 100 -n 10 -M 16000 -aS 0.8 \
            -T ${task.cpus} \
            -o ${prefix}.clustered_80_aS0.8
        
        clstr2txt.pl \
            ${prefix}.clustered_80_aS0.8.clstr \
            > ${prefix}.clustered_80_aS0.8.clstr.txt
        
        mv ${prefix}.clustered_80_aS0.8 ${prefix}.clustered_80_aS0.8.fa
        mv ${prefix}.clustered* $outdir

    """

}
