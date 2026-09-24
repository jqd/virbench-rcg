// MMSEQS2 annotation process doesn't require cdhit input, just the abundance but we include it for consistency
process MMSEQS2_ANNOTATE {
    conda "${params.unified_conda_env}"

    tag "mmseqs2 $sample_id"
    label = 'io_mem'

    publishDir "${params.outdir}/8_mmseqs2/${tool}/${sample_id}", mode: 'copy', overwrite: true, pattern: "${outdir}/${sample_id}.${tool}*"

    publishDir "${params.outdir}/8_mmseqs2/${tool}/log_files/${tool}/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs), path(abundance), path(cdhit_clusters), val(tool)

    output:
    tuple val(sample_id), path("${outdir}/${sample_id}.${tool}.mmseqsTaxaRst_abundance.tsv"), path(abundance), path(cdhit_clusters), val(tool), emit: results
    path("${outdir}/*.tsv")
    path("${outdir}/${sample_id}.${tool}.mmseqsTaxaRst_report.txt")
    path("${outdir}/${sample_id}.${tool}.mmseqsTaxaRst_report.html")
    path '.command.*'

    script:
    outdir = "${tool}_${sample_id}"
    """
    # module load mmseqs2/15-6f452
    mkdir $outdir
    
    # Copy abundance file to output directory
    cp ${abundance} ./$outdir/
    mmseqs createdb $contigs ./$outdir/qry
    #mmseqs createindex ./$outdir/qry tmp

    mmseqs taxonomy ./$outdir/qry ${params.mmseqs2_db} ./$outdir/mmseqsTaxaRst tmp \
     --tax-lineage 1 \
     --majority 0.4 \
     --vote-mode 1 \
     --lca-mode 3 \
     --orf-filter 0 \
     --split-memory-limit 160G \
     --threads 8

    # report
    mmseqs createtsv ./$outdir/qry ./$outdir/mmseqsTaxaRst ./$outdir/${sample_id}.${tool}.mmseqsTaxaRst.tsv
    mmseqs taxonomyreport ${params.mmseqs2_db} ./$outdir/mmseqsTaxaRst ./$outdir/${sample_id}.${tool}.mmseqsTaxaRst_report.txt --report-mode 0
    mmseqs taxonomyreport ${params.mmseqs2_db} ./$outdir/mmseqsTaxaRst ./$outdir/${sample_id}.${tool}.mmseqsTaxaRst_report.html --report-mode 1

    # Parse mmseqs2 output and merge with abundance
    python3 ${projectDir}/bin/mmseqs2_merge_abundance.py \
        --mmseqs2 ./$outdir/${sample_id}.${tool}.mmseqsTaxaRst.tsv \
        --abundance ${abundance} \
        --output ./$outdir/${sample_id}.${tool}.mmseqsTaxaRst_abundance.tsv

    """
}
// less 23_viruses_test.megahit.mmseqsTaxaRst.tsv | cut -f1-2 | less
// k141_139        10359
// k141_47 0
// k141_256        122929
// k141_277        10359
// k141_298        10359