process BLAST_ANNOTATE_COMBINED {
    conda "${params.unified_conda_env}"
    tag "$sample_id"
    label = 'mem_medium'
    
    publishDir "${params.outdir}/6_blast/${tool}/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}.${tool}*.tsv"

    publishDir "${params.outdir}/6_blast/${tool}/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs), path(abundance), path(cdhit_clusters), val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.${tool}_blast_top_hit.tsv"), val(tool), emit: simple_results
    tuple val(sample_id), path("${sample_id}.${tool}_blast_with_abundance.tsv"), val(tool), emit: processed_results
    path("${sample_id}.${tool}.blast.tsv")
    path("${sample_id}.${tool}_blast_processed.tsv")
    path '.command.*'

    script:
    """
    # module load blast/2.9.0
    
    # Run BLAST
    blastn \
     -db  ${params.blast_db} \
     -query $contigs \
     -out ${sample_id}.${tool}.blast.tsv \
     -evalue ${params.blast_evalue} \
     -num_threads ${task.cpus} \
     -outfmt '7 std sskingdom staxids sscinames stitle'

    # Process 1: Get top hit per contig (first hit in BLAST output) - Simple annotation
    grep -v '^#' ${sample_id}.${tool}.blast.tsv | \
        awk 'BEGIN{OFS="\\t"; print "contig_id","taxid"} 
        !seen[\$1]++ {
            # Remove leading "0;" and split by ";" to get first non-zero taxid
            taxid = \$14
            sub(/^0;/, "", taxid)
            split(taxid, arr, ";")
            print \$1, arr[1]
        }' \
        > ${sample_id}.${tool}_blast_top_hit.tsv

    # Process 2: Full processing with Python script
    python ${projectDir}/bin/blast_parser_outfmt7.py \
        --input ${sample_id}.${tool}.blast.tsv \
        --evalue ${params.blast_evalue} \
        --number_bitscores ${params.blast_number_bitscores} \
        ${params.taxonomy ? "--taxonomy ${params.taxonomy}" : "--ncbi_taxonomy_dir ${params.ncbi_taxonomy_dir}"} \
        --output_prefix ${sample_id}.${tool}_blast_processed

    # Process 3: Merge with abundance and cdhit clusters
    python ${projectDir}/bin/blast_merge_abundance.py \
        --blast ${sample_id}.${tool}_blast_processed.tsv \
        --abundance ${abundance} \
        --cdhit ${cdhit_clusters} \
        --output ${sample_id}.${tool}_blast_with_abundance.tsv
    
    """
}
