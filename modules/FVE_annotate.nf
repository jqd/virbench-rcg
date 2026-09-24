process FVE_ANNOTATE {
  conda "${params.unified_conda_env}"
    tag "FVE on $sample_id"
    label = 'io_mem'

    publishDir "${params.outdir}/9_FVE/${sample_id}", 
      mode: 'copy', 
      overwrite: true, 
      pattern: "${outdir}/*.tsv", 
      saveAs: { filename -> filename.split('/').last() }

    publishDir "${params.outdir}/9_FVE/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2) 

    output:
    tuple val(sample_id), path("${outdir}/${tool}_${sample_id}.abundance.tsv"), emit: results
    path("${outdir}/FastViromeExplorer-final-sorted-abundance.tsv")
    path '.command.*'


    script:
    tool = "FVE"
    outdir = "${tool}_${sample_id}"
    """
    #module load java/1.8.0_181
    mkdir $outdir
   

    java -cp "${params.FVE_dir}/bin" FastViromeExplorer -1 $fq1 -2 $fq2 -i ${params.FVE_db} -l ${params.FVE_vlist} -o $outdir
    echo -e "Taxon\tAbundance" > $outdir/${tool}_${sample_id}.abundance.tsv
    sed '1d' $outdir/FastViromeExplorer-final-sorted-abundance.tsv | \
      awk -F'\t' '{split(\$1,a,"|"); sum[a[1]]+=\$4} END{for(i in sum) printf "%s\\t%.4f\\n", i, sum[i]}' | \
      sed 's/taxon://' | sort -k2,2 -nr >> $outdir/${tool}_${sample_id}.abundance.tsv
    #rm $outdir/FastViromeExplorer-reads-mapped-sorted*.bam
    #rm $outdir/FastViromeExplorer*.sam	
    """
}

