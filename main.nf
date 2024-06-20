#!/usr/bin/env nextflow

/*************************************************
* DSL2
**************************************************/
//enable DSL2 syntax
nextflow.enable.dsl=2

/*************************************************
* VERSION
**************************************************/
version = '0.0.9000'

/*************************************************
* FUNCTIONS
**************************************************/

/*
 Usage:
    nextflow run main.nf --input <dir> --transcriptome <file> --salmon_index <dir>
*/
def helpMessage() {
  log.info """
        
        =========================================
        nf-rnaseq: a nextflow rna-seq pipeline
        =========================================
        version ${version}
        
        Usage:
        The typical command for running the pipeline is as follows:

            nextflow run nf-rnaseq -profile <conda/docker/singularity> --modules salmon --input /path/to/samples --transcriptome /path/to/transcriptome --salmon_index /path/to/index [OPTIONS]

        Mandatory arguments:
            --input             DIRPATH             Folder containing FASTQ files, or file with
                                                    4 columns: id, read_1, read_2, library_type
            -profile            STRING              Container to use, one of conda, singularity, or docker.


        Optional arguments:
            --modules           STRING              The pipeline modules to run. By default, all available modules will run 
                                                    ('--modules fastqc,cutadapt,fastqc_cutadapt,salmon,star,multiqc').
            --filext            STRING              Extension of input files (default: fq.gz)
            --suffix1           STRING              Suffix of first file in paired reads (default: _1)
            --suffix2           STRING              Suffix of second file in paired reads (default: _2)
            --concatenate       BOOLEAN             Whether to concatenate input files when multiple files 
                                                    per sample id are found (e.g., files from different lanes)
            --prefix            STRING              Regular expression used to identify groups of multiple
                                                    files to concatenate (e.g., --prefix LANE(\d+)_)
            --species           STRING              Species of the samples (e.g., --species hsapiens). 
                                                    This parameter is used to create the output sub-folders 
                                                    and to download genome/transcriptome data (if required)
            --refdir            DIRPATH             Folder with reference transcriptome and (optional) genome
            --decoys           [FILEPATH]           File containing a set of decoy sequences. If the parameter is
                                                    provided without value (i.e. --decoys), a set of decoys
                                                    will be computed from the transcriptome and genome files.
            --genome            FILEPATH            FASTA file containing the genome (can be a gzip file)
            --transcriptome     FILEPATH            FASTA file containing the transcriptome (can be a gzip file)
            --gtf               FILEPATH            Gene Transfer Format file (can be a gzip file)
            --genemap          [FILEPATH]           File containing a mapping of transcripts to genes used if the salmon 
                                                    module is selected.  If the parameter is provided without a
                                                    value (i.e. --genemap), and a GTF file is provided in input,
                                                    mapping is attempted.
            --salmon_index      [DIRPATH]           Folder containing the salmon index. If the parameter is provided
                                                    without a value (i.e. --salmon_index), an index will be generated
                                                    from the transcriptome and genome files. This is a mandatory 
                                                    argument if salmon is selected.
            --star_index        [DIRPATH]           Folder containing the STAR index. If the parameter is provided
                                                    without a value (i.e. --star_index), an index will be generated
                                                    from the genome and GTF files. This is a mandatory argument
                                                    if star is selected.                                            
            --salmon_libtype    STRING              Library type, used for salmon quantification (default: 'A')
            --stranded          STRING              A single integer denoting the strandedness of the reads for 
                                                    featureCounts quantification. Possible values include:
                                                    0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                                                    This is a mandatory argument if star alignment is selected.
            --multiqc_config    FILEPATH            Config yaml file for MultiQC
            --outdir            DIRPATH             Output directory (default: ./results)
            --cachedir          DIRPATH             Provide a centralised cache directory for containers (default: ./work)
            --verbose                               Whether to report extra information on progress
            --help                                  Print this usage statement
            --max_cpus          STRING              Maximum amount of allowed cpus (default: 7)
            --max_memory        STRING              Maximum amount of allowed memory (default: '30.GB')
            --max_time          STRING              Maximum amount of execution time (default: '48.h')
        """.stripIndent()
}

def initialLogMessage() {
    log.info """

        =========================================
        nf-rnaseq: a nextflow rna-seq pipeline
        =========================================
        version ${version}
        
        Usage:
        nextflow run main.nf --input <PATH> --refdir <PATH> --salmon_index <PATH>
    """.stripIndent()

    log.info """

    -----------------------------------------
    CURRENT SETTINGS
    -----------------------------------------
    """.stripIndent()

    log.info "\nCONFIGURATION ---------------------------"
    if (params.modules       ) log.info "modules        = ${params.modules}"
    if (params.multiqc_config) log.info "multiqc_config = ${params.multiqc_config}"
    if (params.verbose       ) log.info "verbose        = ${params.verbose}"
    if (params.cachedir      ) log.info "cachedir       = ${params.cachedir}"

    log.info "\nINPUT/OUTPUT ----------------------------"
    if (params.input         ) log.info "input          = ${params.input}"
    if (params.prefix        ) log.info "prefix         = ${params.prefix}"
    if (params.suffix1       ) log.info "suffix1        = ${params.suffix1}"
    if (params.suffix2       ) log.info "suffix2        = ${params.suffix2}"
    if (params.filext        ) log.info "filext         = ${params.filext}"
    if (params.outdir        ) log.info "outdir         = ${params.outdir}"

    log.info "\nREFERENCES ------------------------------"   
    if (params.species       ) log.info "species        = ${params.species}"
    if (params.genome        ) log.info "genome         = ${params.genome}"
    if (params.transcriptome ) log.info "transcriptome  = ${params.transcriptome}"
    if (params.gtf           ) log.info "GTF            = ${params.gtf}"
    if (params.decoys        ) log.info "decoys         = ${params.decoys}"
    if (params.genemap       ) log.info "genemap        = ${params.genemap}"
    if (params.salmon_index  ) log.info "salmon_index   = ${params.salmon_index}"
    if (params.star_index    ) log.info "star_index     = ${params.star_index}"

    log.info "\nQUANTIFICATION---------------------------"   
    if (params.salmon_libtype) log.info "salmon_libtype = ${params.salmon_libtype}"
    if (params.stranded      ) log.info "strandedness   = ${params.stranded}"

    log.info "\nRESOURCES -------------------------------"   
    if (params.max_memory    ) log.info "max_memory     = ${params.max_memory}"
    if (params.max_cpus      ) log.info "max_cpus       = ${params.max_cpus}"
    if (params.max_time      ) log.info "max_time       = ${params.max_time}"
        
    log.info """

    -----------------------------------------------
    END SETTINGS
    -----------------------------------------------
    """.stripIndent()

    // log.info """
    //     /---------------------------------------------\\
    //     |
    //     | input             = ${params.input}
    //     | prefix            = ${params.prefix}
    //     | suffix1           = ${params.suffix1}
    //     | suffix2           = ${params.suffix2}
    //     | filext            = ${params.filext}
    //     | outdir            = ${params.outdir}
    //     | multiqc_config    = ${params.multiqc_config}
    //     |
    //     | species           = ${params.species}
    //     | genome            = ${params.genome}
    //     | transcriptome     = ${params.transcriptome}
    //     | gtf               = ${params.gtf}
    //     | genemap           = ${params.genemap}
    //     | salmon_index      = ${params.salmon_index}
    //     | star_index        = ${params.star_index}
    //     \\---------------------------------------------/

    // """.stripIndent()
}


/*************************************************
* CONFIG FILES
**************************************************/
//MultiQC default config file
def multiqc_config_fpath = "$projectDir/config/multiqc_config.yaml"

/*************************************************
* PIPELINE PARAMETERS
**************************************************/

// Get the modules to run as elements of array
def modules_to_run = params.modules ? "${params.modules}".split(',') : []

/*************************************************
* IMPORT MODULES
**************************************************/
include { FASTQC                    } from './modules/fastqc'
include { MULTIQC                   } from './modules/multiqc'
include { CUTADAPT	            } from './modules/cutadapt'
include { FASTQC as FASTQC_CUTADAPT } from './modules/fastqc'
include { SALMON_QUANT              } from './modules/salmon/quant'
include { TXIMPORT_SALMON           } from './modules/tximport/salmon'
include { STAR_QUANT                } from './modules/star/quant'
include { FEATURE_COUNTS            } from './modules/featureCounts'
include { SUMMARISE_FEATURECOUNTS   } from './modules/summarise_featurecounts'
include { SOFTWARE_VERSIONS         } from './modules/custom/swversions'

/*************************************************
* IMPORT SUBWORKFLOWS
**************************************************/
include { CHECK_AND_PREPARE_INPUT } from './subworkflows/check_and_prepare_input'
include { CONCATENATE_READS       } from './subworkflows/concatenate_reads' 

/*************************************************
* RUN WORKFLOW
**************************************************/
workflow {
    
    //---------------------------------------------
    // LOG
    //---------------------------------------------
    if(params.help) {
        helpMessage()
        exit 0
    } else {
        if(params.verbose) initialLogMessage()
    }

    //---------------------------------------------
    // SOFTWARE VERSIONS
    //---------------------------------------------
    ch_versions = Channel.empty()

    //---------------------------------------------
    // CHECK INPUT
    //---------------------------------------------
    // Check input is provided
    CHECK_AND_PREPARE_INPUT()

    // Set channels
    ch_reads         = CHECK_AND_PREPARE_INPUT.out.reads
    ch_salmon_index  = CHECK_AND_PREPARE_INPUT.out.salmon_index
    ch_star_index    = CHECK_AND_PREPARE_INPUT.out.star_index
    ch_transcriptome = CHECK_AND_PREPARE_INPUT.out.transcriptome
    ch_gtf           = CHECK_AND_PREPARE_INPUT.out.gtf
    ch_genemap       = CHECK_AND_PREPARE_INPUT.out.genemap
    
    if(params.verbose) ch_reads.view(it -> "[READS] $it")
    if(params.verbose) ch_salmon_index.view(it -> "[INDEX] $it")
    if(params.verbose) ch_star_index.view(it -> "[INDEX] $it")
    if(params.verbose) ch_transcriptome.view(it -> "[TRANSCRIPTOME] $it")
    if(params.verbose) ch_genemap.view(it -> "[GENEMAP] $it")
    
    // Update versions
    ch_versions = ch_versions.mix(CHECK_AND_PREPARE_INPUT.out.versions)

    //---------------------------------------------
    // FastQC
    //---------------------------------------------
    
    if(modules_to_run.contains('fastqc')){
        ch_fastqc_dir = 'fastqc'
        FASTQC(
            ch_reads,
            ch_fastqc_dir
        )
        // Update versions
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    if(params.concatenate){
        CONCATENATE_READS(
            ch_reads
        )
        // Update
        ch_reads = CONCATENATE_READS.out.reads
        ch_versions = ch_versions.mix(CONCATENATE_READS.out.versions)
    }

    //---------------------------------------------
    // Adapter trimming
    //---------------------------------------------    
    if(modules_to_run.contains('cutadapt')){
	CUTADAPT(
            ch_reads
        )
        // Update 
	ch_reads = CUTADAPT.out.reads
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    if(modules_to_run.contains('fastqc_cutadapt')){
        if(modules_to_run.contains('fastqc') && !modules_to_run.contains('cutadapt')){
            exit 1, 'ERROR: The fastqc and fastqc_cutadapt modules may not be specified together unless the cutadapt module is also present'
        }

        if(modules_to_run.contains('cutadapt') | params.input.endsWith('cutadapt/')){
            ch_fastqc_dir = 'fastqc_cutadapt'
            FASTQC_CUTADAPT(
                ch_reads,
                ch_fastqc_dir
            )
            // Update versions
            ch_versions = ch_versions.mix(FASTQC_CUTADAPT.out.versions.first())
        } else {
            exit 1, 'ERROR: The fastqc_cutadapt module is specified without the cutadapt module, and the input files are not stored in a cutadapt/ directory. Update the input directory and the suffix parameters to feed the trimmed files into fastqc_cutadapt directly.'
        }
    }

    //---------------------------------------------
    // Pseudo-alignment and quantification
    //---------------------------------------------
    if(modules_to_run.contains('salmon')){
        // LOG
        // ch_salmon_index.ifEmpty([]).view(it -> "[SALMON INDEX] $it")
        // ch_transcriptome.ifEmpty([]).view(it -> "[TRANSCRIPTOME] $it")
        // ch_genemap.ifEmpty([]).view(it -> "[GENEMAP] $it")
        // Channel.value('A').view()
        SALMON_QUANT(
            ch_reads,
            ch_salmon_index.collect(),
            ch_transcriptome.collect(),
            ch_genemap.collect().ifEmpty(file('NO_FILE')),
            Channel.value(params.salmon_libtype)
        )
        // Update channels
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

        // Summarisation at gene-level
        if(ch_genemap){
            TXIMPORT_SALMON(
                SALMON_QUANT.out.quant.collect{it[1]},
                ch_genemap.collect()
            )
            // Update channels
            ch_versions = ch_versions.mix(TXIMPORT_SALMON.out.versions)
        }

    }

    //---------------------------------------------
    // STAR-alignment and quantification
    //---------------------------------------------
    
    if(modules_to_run.contains('star')){
        STAR_QUANT(
            ch_reads,
            ch_star_index.collect(),
            ch_gtf.collect()
        )
        // Update
        ch_bam = STAR_QUANT.out.bam_channel
        ch_versions = ch_versions.mix(STAR_QUANT.out.versions.first())

        FEATURE_COUNTS(
            ch_bam,
            ch_gtf.collect(),
            Channel.value(params.stranded)
        )

        // Update
        ch_featurecounts = FEATURE_COUNTS.out.featurecounts
        ch_versions = ch_versions.mix(FEATURE_COUNTS.out.versions.first())

        SUMMARISE_FEATURECOUNTS(
            ch_featurecounts.collect{it[1]}
        )

        // Update
        ch_versions = ch_versions.mix(SUMMARISE_FEATURECOUNTS.out.versions)
    }


    //---------------------------------------------
    // Software versions
    //---------------------------------------------
    SOFTWARE_VERSIONS(
        ch_versions.unique().collectFile(name: 'merged_software_versions.yml')
    )


    //---------------------------------------------
    // MultiQC
    //---------------------------------------------
    if(modules_to_run.contains('multiqc')){
        //Channels
        ch_multiqc_config = Channel.fromPath( multiqc_config_fpath, checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
        //Choose one
        ch_multiqc_config = params.multiqc_config ? ch_multiqc_custom_config : ch_multiqc_config

        //Run module
        MULTIQC(
            ch_multiqc_config,
            SOFTWARE_VERSIONS.out.multiqc_yml.collect().ifEmpty(file('NO_FILE'))
        )
        // Update channels
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

}

//---------------------------------------------
// COMPLETION
//---------------------------------------------
workflow.onComplete {
    def msg = """\
        
        
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
    print msg
}
