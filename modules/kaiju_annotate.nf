// Kaiju annotation process doesn't require cdhit input, just the abundance but we include it for consistency
process KAIJU_ANNOTATE {
    conda "${params.unified_conda_env}"
    tag "kaiju on $sample_id"
    label = 'io_mem'

    publishDir "${params.outdir}/7_kaiju/${tool}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.tsv'

    publishDir "${params.outdir}/7_kaiju/${tool}/log_files/${tool}/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs), path(abundance), path(cdhit_clusters), val(tool)

    output:
    tuple val(sample_id), path("${sample_id}.${tool}.kaiju_abundance.tsv"), path(abundance), path(cdhit_clusters), val(tool), emit: results
    path("${sample_id}.${tool}.kaiju.tsv")
    path '.command.*'

    script:
    """
    # Run Kaiju on contigs (requires Kaiju and database)
    kaiju -z 5 -v \
          -t ${params.kaiju_tax_path}/nodes.dmp \
          -f ${params.kaiju_db} \
          -i $contigs \
          -a greedy \
          -o ${sample_id}.${tool}.kaiju.tsv

    # Parse kaiju output and merge with abundance
    python3 ${projectDir}/bin/kaiju_merge_abundance.py \
        --kaiju ${sample_id}.${tool}.kaiju.tsv \
        --abundance ${abundance} \
        --output ${sample_id}.${tool}.kaiju_abundance.tsv
    """
}

// Example kaiju output:
// less 23_viruses_test.megahit.kaiju.tsv | cut -f1-3
// C       k141_235        10298
// C       k141_236        10359
// C       k141_119        10359
// C       k141_120        10298
// C       k141_121        10298
// C       k141_122        10359

//map_back_to_contigs output:
// output:
//   tuple val(sample_id), path(contigs), path("${sample_id}_${tool}.abundance.tsv"), val(tool), emit: abundance
// Example abundance output:
// k141_709        92
// k141_1  6
// k141_354        16
// k141_1063       22
// k141_1065       1
// k141_2  76

//cdhit output:
// output:
//   tuple val(sample_id), path("${outdir}/${prefix}.clustered*.txt"), val(tool), emit: clusters

// EXAMPLE cdhit clustered output:
// id      clstr   clstr_size      length  clstr_rep       clstr_iden      clstr_cov
// k141_1346       0       5       18131   1       100     100%
// k141_739        0       5       3407    0       87.03%  18%
// k141_758        0       5       1466    0       81.58%  8%
// k141_19 0       5       1002    0       87.13%  5%
