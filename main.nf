
nextflow.enable.dsl=2

include { FASTP } from './modules/fastp.nf'
include { REMOVE_HOST_READS } from './modules/remove_host_reads.nf'
include { ASSEMBLE_MEGAHIT } from './modules/assembly_megahit.nf'
include { ASSEMBLE_METASPADES } from './modules/assembly_metaspades.nf'

include { MAP_BACK_TO_CONTIGS as MAP_BACK_TO_CONTIGS_1 } from './modules/map_back_to_contigs.nf'
include { MAP_BACK_TO_CONTIGS as MAP_BACK_TO_CONTIGS_2 } from './modules/map_back_to_contigs.nf'

include { CDHIT as CDHIT_CLUSTER_1 } from './modules/cdhit.nf'
include { CDHIT as CDHIT_CLUSTER_2 } from './modules/cdhit.nf'

include { BLAST_ANNOTATE_COMBINED as BLAST_ANNOTATE_1 } from './modules/blast_annotate_combined.nf'
include { BLAST_ANNOTATE_COMBINED as BLAST_ANNOTATE_2 } from './modules/blast_annotate_combined.nf'

include { KAIJU_ANNOTATE as KAIJU_ANNOTATE_1 } from './modules/kaiju_annotate.nf'
include { KAIJU_ANNOTATE as KAIJU_ANNOTATE_2 } from './modules/kaiju_annotate.nf'

include { MMSEQS2_ANNOTATE as MMSEQS2_ANNOTATE_1 } from './modules/mmseqs2_annotate.nf'
include { MMSEQS2_ANNOTATE as MMSEQS2_ANNOTATE_2 } from './modules/mmseqs2_annotate.nf'

include { FVE_ANNOTATE } from './modules/FVE_annotate.nf'
include { VIROMESCAN_ANNOTATE } from './modules/viromescan_annotate.nf'

include { KRAKEN2_ANNOTATE } from './modules/kraken2_annotate.nf'


include { CALCULATE_METRICS as CALCULATE_METRICS_1 } from './modules/calculate_metrics.nf'
include { CALCULATE_METRICS as CALCULATE_METRICS_2 } from './modules/calculate_metrics.nf'
include { GATHER_RESULTS as GATHER_RESULTS_1 } from './modules/gather_results.nf'
include { GATHER_RESULTS as GATHER_RESULTS_2 } from './modules/gather_results.nf'

params.help = params.help ?: false

def printHelp() {
    log.info """
    DESCRIPTION:
    virbench-rcg workflow for multi-tool virome profiling and comparative benchmarking
    across database and assembly settings.

    Raw sequencing reads are first preprocessed by quality control and host-read removal.
    The workflow then branches into read-based profiling (Kraken2, ViromeScan, and
    FastViromeExplorer) and contig-based profiling, where reads are assembled (MEGAHIT
    or metaSPAdes), clean reads are mapped back to contigs for abundance estimation,
    and contigs are annotated (BLAST, Kaiju, and MMseqs2 taxonomy).

    Outputs from both branches are evaluated under two database configurations (uniform
    RefSeq-CG database and tool-specific databases). Final benchmarking compares tool 
    performance using precision, recall, F1 score, and L2 distance.

    Usage:
        nextflow -C nextflow.config run main.nf [options]

    Core options:
        --samples <file>                 Sample sheet TSV (default: samples.tsv)
        --outdir <dir>                   Output directory (default: results)
        --assembler <name>               megahit | metaspades | both (default: both)
        --skip_fastp <true|false>        Skip read quality trimming and use input reads directly (default: false)
        --skip_metrics <true|false>      Skip metrics calculation (default: false)
        --skip_remove_host <true|false>  Skip host-read removal and use fastp reads directly (default: false)

    Read preprocessing:
        --host <bowtie2_index_prefix>    Host genome Bowtie2 index prefix

    Annotation database options:
        --blast_db <path>                BLAST nucleotide database
        --blast_evalue <value>           BLAST e-value threshold (default: 1e-5)
        --blast_number_bitscores <int>   Number of BLAST bitscore hits kept (default: 5)
        --kaiju_db <path>                Kaiju database (.fmi)
        --kaiju_tax_path <path>          Kaiju taxonomy path
        --mmseqs2_db <path>              MMseqs2 database
        --kraken2_db <path>              Kraken2 database
        --FVE_db <path>                  FastViromeExplorer index
        --FVE_vlist <path>               FastViromeExplorer virus list
        --FVE_dir <path>                 FastViromeExplorer installation directory
        --FVE_default_db <true|false>    FVE uses default DB accession-to-taxon mode
        --FVE_acc2taxon <path>           FVE accession-to-taxon mapping file
        --viromescan_db <name|path>      ViromeScan database identifier/path
        --viromescan_default_db <true|false> ViromeScan uses default DB accession-to-taxon mode
        --viromescan_acc2taxon <path>    ViromeScan accession-to-taxon mapping file

    CD-HIT options:
        --cdhit_identity <float>         Sequence identity threshold (default: 0.95)
        --cdhit_word_length <int>        CD-HIT word length (default: 10)

    Taxonomy / metrics references:
        --ncbi_taxonomy_dir <path>       NCBI taxonomy directory for lineage parsing
        --taxonomy <path>                Optional precomputed taxonomy lineage file
        --ground_truth <path>            Ground-truth abundance/lineage table

    Profiles:
        -profile standard                Run locally (executor profile)
        -profile lsf                     Run on LSF cluster (executor profile)
        -profile db_default              Use tool-specific default databases
        -profile db_refseq               Use uniform RefSeq-CG databases

    HPC profile notes:
        Select the scheduler profile that matches your HPC system. Review
        queue names and scheduler settings in nextflow.config and
        profiles/schedulers.config before running at your institution.

    Other:
        --help                           Show this help message and exit

    Example:
        # Default tool-specific DB profile
        nextflow -C nextflow.config run main.nf \
            --samples samples.tsv --assembler both --skip_metrics false \
            -profile lsf,db_default -resume

        # Uniform RefSeq-CG DB profile
        nextflow -C nextflow.config run main.nf \
            --samples samples.tsv --assembler both --skip_metrics false \
            -profile lsf,db_refseq -resume
    """.stripIndent()
}

if (params.help) {
        printHelp()
        System.exit(0)
}

def validateDatabaseParams() {
    def requiredCoreParams = [
        'outdir',
        'host'
    ]
    if (!params.skip_metrics) {
        requiredCoreParams += ['ncbi_taxonomy_dir', 'ground_truth']
    }

    def requiredDbParams = [
        'blast_db',
        'kaiju_db',
        'kaiju_tax_path',
        'mmseqs2_db',
        'FVE_db',
        'FVE_vlist',
        'viromescan_db',
        'kraken2_db'
    ]

    def missing = []
    missing.addAll(requiredCoreParams.findAll { !params[it] })
    missing.addAll(requiredDbParams.findAll { !params[it] })

    if (params.FVE_default_db == null)
        missing << 'FVE_default_db'
    if (params.viromescan_default_db == null)
        missing << 'viromescan_default_db'

    if ("${params.FVE_default_db}".toLowerCase() == 'true' && !params.FVE_acc2taxon)
        missing << 'FVE_acc2taxon'
    if ("${params.viromescan_default_db}".toLowerCase() == 'true' && !params.viromescan_acc2taxon)
        missing << 'viromescan_acc2taxon'

    missing = missing.unique()
    if (missing) {
        error """
        Missing required parameters: ${missing.join(', ')}
        Select a DB profile, for example:
                    -profile lsf,db_default
                    -profile lsf,db_refseq
        Set profile paths by environment variables or pass values explicitly with --<param>.
        Or provide values explicitly with --<param> options.
        """.stripIndent().trim()
    }
}

validateDatabaseParams()

workflow {

    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.fq1), file(row.fq2)) }

    fastp_out = params.skip_fastp
        ? [reads: samples_ch]
        : FASTP(samples_ch)

    clean_reads = (params.skip_remove_host)
        ? [reads: fastp_out.reads]
        : REMOVE_HOST_READS(fastp_out.reads)

// BLAST annotation workflow has two approaches:
// 1. Extract top hit only - no integration CD-HIT clustering (optional)
// 2. Full processing - use custom Python script to parse BLAST output and merge with CD-HIT results (default)

    // Run selected assembler(s)
    if (params.assembler == 'megahit' || params.assembler == 'both') {
        megahit_out = ASSEMBLE_MEGAHIT(clean_reads.reads)
        merged_megahit_ch = clean_reads.reads.combine(megahit_out.contigs, by: 0).map { it + ['megahit'] }
        megahit_abun_ch = MAP_BACK_TO_CONTIGS_1(merged_megahit_ch)
        megahit_cdhit_ch = CDHIT_CLUSTER_1(megahit_out.contigs.map { sample_id, contigs -> tuple(sample_id, contigs, 'megahit') })
        
        // Prepare input for annotation tools (contigs + abundance + cdhit)
        megahit_annotate_input = megahit_abun_ch.abundance
            .map { sample_id, contigs, abundance, tool -> tuple(sample_id, contigs, abundance, tool) }
            .combine(megahit_cdhit_ch.clusters.map { sample_id, cdhit_file, tool -> tuple(sample_id, cdhit_file) }, by: 0)
            .map { sample_id, contigs, abundance, tool, cdhit_file -> tuple(sample_id, contigs, abundance, cdhit_file, tool) }
        
        // Run annotation tools
        blast_megahit_ch = BLAST_ANNOTATE_1(megahit_annotate_input)
        kaiju_megahit_ch = KAIJU_ANNOTATE_1(megahit_annotate_input)
        mmseqs2_megahit_ch = MMSEQS2_ANNOTATE_1(megahit_annotate_input)
    }
    
    if (params.assembler == 'metaspades' || params.assembler == 'both') {
        metaspades_out = ASSEMBLE_METASPADES(clean_reads.reads)
        merged_metaspades_ch = clean_reads.reads.combine(metaspades_out.contigs, by: 0).map { it + ['metaspades'] }
        metaspades_abun_ch = MAP_BACK_TO_CONTIGS_2(merged_metaspades_ch)
        metaspades_cdhit_ch = CDHIT_CLUSTER_2(metaspades_out.contigs.map { sample_id, contigs -> tuple(sample_id, contigs, 'metaspades') })
        
        // Prepare input for annotation tools (contigs + abundance + cdhit)
        metaspades_annotate_input = metaspades_abun_ch.abundance
            .map { sample_id, contigs, abundance, tool -> tuple(sample_id, contigs, abundance, tool) }
            .combine(metaspades_cdhit_ch.clusters.map { sample_id, cdhit_file, tool -> tuple(sample_id, cdhit_file) }, by: 0)
            .map { sample_id, contigs, abundance, tool, cdhit_file -> tuple(sample_id, contigs, abundance, cdhit_file, tool) }
        
        // Run annotation tools
        blast_metaspades_ch = BLAST_ANNOTATE_2(metaspades_annotate_input)
        kaiju_metaspades_ch = KAIJU_ANNOTATE_2(metaspades_annotate_input)
        mmseqs2_metaspades_ch = MMSEQS2_ANNOTATE_2(metaspades_annotate_input)
    }

    FVE_ch = FVE_ANNOTATE(clean_reads.reads)
    VIROMESCAN_ch = VIROMESCAN_ANNOTATE(clean_reads.reads)
    KRAKEN2_ch = KRAKEN2_ANNOTATE(clean_reads.reads)

    // Prepare channels for metrics calculation based on assembler selection
    if (params.assembler == 'megahit' || params.assembler == 'both') {
        // Combine all tool results by sample_id for megahit
        megahit_metrics_input = blast_megahit_ch.processed_results
            .map { sample_id, file, tool -> tuple(sample_id, file, tool) }
            .join(kaiju_megahit_ch.results.map { sample_id, file, abundance, cdhit, tool -> tuple(sample_id, file) })
            .join(mmseqs2_megahit_ch.results.map { sample_id, file, abundance, cdhit, tool -> tuple(sample_id, file) })
            .join(KRAKEN2_ch.results)
            .join(FVE_ch.results)
            .join(VIROMESCAN_ch.results)

        // Gather all tool results into a single folder per sample
        megahit_gathered_results = GATHER_RESULTS_1(megahit_metrics_input)

        // Calculate metrics if not skipped
        if (!params.skip_metrics) {
            megahit_metrics_input_with_refs = megahit_metrics_input
                .map { sample_id, blast, blast_tool, kaiju, mmseqs2, kraken2, fve, viromescan ->
                    tuple(sample_id, blast, blast_tool, kaiju, mmseqs2, kraken2, fve, viromescan, 
                          file(params.ground_truth)) }
            megahit_metrics_ch = CALCULATE_METRICS_1(megahit_metrics_input_with_refs)
        }
    }
    
    if (params.assembler == 'metaspades' || params.assembler == 'both') {
        // Combine all tool results by sample_id for metaspades
        metaspades_metrics_input = blast_metaspades_ch.processed_results
            .map { sample_id, file, tool -> tuple(sample_id, file, tool) }
            .join(kaiju_metaspades_ch.results.map { sample_id, file, abundance, cdhit, tool -> tuple(sample_id, file) })
            .join(mmseqs2_metaspades_ch.results.map { sample_id, file, abundance, cdhit, tool -> tuple(sample_id, file) })
            .join(KRAKEN2_ch.results)
            .join(FVE_ch.results)
            .join(VIROMESCAN_ch.results)

        // Gather all tool results into a single folder per sample
        metaspades_gathered_results = GATHER_RESULTS_2(metaspades_metrics_input)

        // Calculate metrics if not skipped
        if (!params.skip_metrics) {
            metaspades_metrics_input_with_refs = metaspades_metrics_input
                .map { sample_id, blast, blast_tool, kaiju, mmseqs2, kraken2, fve, viromescan ->
                    tuple(sample_id, blast, blast_tool, kaiju, mmseqs2, kraken2, fve, viromescan, 
                          file(params.ground_truth)) }
            metaspades_metrics_ch = CALCULATE_METRICS_2(metaspades_metrics_input_with_refs)
        }
    }
    
}
