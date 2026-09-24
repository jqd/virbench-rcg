process VIROMESCAN_ANNOTATE {
  conda "${params.viromescan_conda_env}"
    tag "viromescan on $sample_id"
    label = 'mem_veryhigh'

    publishDir "${params.outdir}/10_viromescan/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}*.tsv"

    publishDir "${params.outdir}/10_viromescan/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2) 

    output:
    tuple val(sample_id), path("${sample_id}.viromescan.species_abundance.tsv"), emit: results
    path("${sample_id}*.tsv")
    path '.command.*'


    script:
    tool = "viromescan"
    outdir = "${tool}_${sample_id}"
    """
    viromescan -p 3 -m ${params.viromescan_conda_env} -d ${params.viromescan_db} -1 $fq1 -2 $fq2 -o $outdir
    echo -e "Taxon\tAbundance" > ${sample_id}.viromescan.species_abundance.tsv
    awk '\$3!=0 {split(\$1,a,"|"); key=(a[1]=="gi" ? a[4] : a[2]); s[key]+=\$3} END{for(i in s) print i "\t" s[i]}' ${outdir}/final.genes.txt | \
      sort -k2,2 -nr >> ${sample_id}.viromescan.species_abundance.tsv

    """
}

