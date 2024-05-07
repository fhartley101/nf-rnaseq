#!/usr/bin/env nextflow

/*************************************************
* DSL2
**************************************************/
//enable DSL2 syntax
nextflow.enable.dsl=2

/*************************************************
* FUNCTIONS
**************************************************/
//Create reads input channel from directory path
def create_channel_from_dir(String dirpath, String filext, String suffix1, String suffix2, String prefix){
    
    if (!(prefix instanceof String)){
        prefix = ''
    }
    
    if(file(dirpath).isDirectory() ){
        if (file(dirpath).listFiles().findAll { it.name ==~ /.*${filext}/ }.size() > 0){
            // print "Input files found: ";
            //Check if has paired reads
            if (file(dirpath).listFiles().findAll { it.name ==~ /.*${suffix2}.${filext}/ }.size() > 0){
                //Paired reads
                println "Input files found: paired library";
                reads = "${dirpath}/*{${suffix1},${suffix2}}.fastq.gz"
                single_end = false
            } else {
                //Single reads
                println "Input files found: single library";
                reads = file(dirpath).listFiles().findAll { it.name ==~ /.*.${filext}/ }
                single_end = true
            }
            //Create input channel
            ch_reads = Channel
                .fromFilePairs(reads)
            // Replace prefix, suffix and file extension
            ch_reads = ch_reads.map { 
                it -> [
                    it[0].replaceAll( ~/(${prefix})|(${suffix1})|(((${suffix1})?.${filext}))/, "" ),//regexp: (_1)?.fq.gz
                    it[1]
                ]
            }
            // Group tuple
            ch_reads = ch_reads.groupTuple()
            // Update map
            if(single_end){
                ch_reads = ch_reads.map {
                    it -> 
                        meta = [
                            id: it[0],
                            single_end: true
                        ]
                        // Flatten all array entries
                        reads1 = it[1].flatten()
                        // Collect the reads into an array (multi-lane case), sort the array
                        reads1 = reads1.collect().sort()
                        // Return tuple
                        return [meta, reads1, []]
                }
            } else {
                ch_reads = ch_reads.map {
                    it -> 
                        // Create metadata array
                        meta = [
                            id: it[0],
                            single_end: false
                        ]
                        // Flatten all array entries
                        reads = it[1].flatten()
                        // Divide the reads, collect them into an array, sort the array. Sorting is needed
                        // to make sure that in case of multi-lane input data for the paired-reads is 
                        // matching. For example, let's assume we have 2 lanes for one paired reads. 
                        // We want to have the same order, like: [LANE1_READ_1, LANE2_READ_1] and [LANE1_READ_2, LANE2_READ_2]
                        reads1 = reads.findAll{ it.name ==~ /.*${suffix1}.${filext}/ }.collect().sort()
                        reads2 = reads.findAll{ it.name ==~ /.*${suffix2}.${filext}/ }.collect().sort()
                        // Return tuple
                        return [meta, reads1, reads2]
                }
            }

        } else {
            println "ERROR: no input data with the given extension found";
            System.exit(1)
        }
    } else {
        println "ERROR: input type not supported";
        System.exit(1)
    }

    return ch_reads
}

//Download genome/transcriptome
def downloadFile(String aUrl, String aPath, String aFileName, Boolean printLog){
    // URL
    url = new URL(aUrl)
    // File name
    if(aFileName == null){
        aFileName = file(aUrl).getName()
    }
    // File path
    aFilePath = "${aPath}/${aFileName}"
    // Define new file
    def aFile = new File( aFilePath )
    // Get file
    if ( aFile.exists() ){
        log.info "File found in pipeline folder: delete the file for a new download"
        log.info "File path: ${aFilePath}"
    } else {
        log.info "URL: ${aUrl}"
        log.info "File path: ${aFilePath}"
        // Check if dir exists
        if( !file(aPath).exists() ){
            if(printLog) log.info "The directory where to store reference data DOES NOT exists and will be created"
            (new File( aPath )).mkdirs()
        }
        if(printLog) log.info "Downloading the file"
        url.withInputStream {
            inputStream -> aFile << inputStream
        }
    }

    return aFilePath
}

def getTranscriptomeUrl(species){
    switch(species){
        case ~/(homo_sapiens|hsapiens)/ : return "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
        case "mus_musculus"             : return "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
        case "mus_musculus_fvbnj"       : return "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus_fvbnj/cdna/Mus_musculus_fvbnj.FVB_NJ_v1.cdna.all.fa.gz"
    }
    // throw new IllegalArgumentException("Provided species (" + species + ") is not currently supported.")
    exit 1, "Provided species (" + species + ") is not currently supported."
}

def getGenomeUrl(species){
    switch(species){
        case ~/(homo_sapiens|hsapiens)/ : return "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        case "mus_musculus"             : return "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        case "mus_musculus_fvbnj"       : return "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus_fvbnj/dna/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa.gz"
    }
    // throw new IllegalArgumentException("Provided species (" + species + ") is not currently supported.")
    exit 1, "Provided species (" + species + ") is not currently supported."
}

def getGtfUrl(species){
    switch(species){
        case ~/(homo_sapiens|hsapiens)/ : return "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"
        case "mus_musculus"             : return "https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"
        case "mus_musculus_fvbnj"       : return "https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus_fvbnj/Mus_musculus_fvbnj.FVB_NJ_v1.109.gtf.gz"
    }
    // throw new IllegalArgumentException("Provided species (" + species + ") is not currently supported.")
    exit 1, "Provided species (" + species + ") is not currently supported."
}

def supported_species = ["homo_sapiens", "hsapiens", "mus_musculus", "mus_musculus_fvbnj"]
def supported_modules = ["fastqc", "quant", "multiqc"]
def raw_reads_filext = ['fastq', 'fq', 'fastq.gz', 'fq.gz']
def aligned_reads_filext    = ['bam', 'sam']
def supported_filext  = raw_reads_filext + aligned_reads_filext

def modules_to_run = params.modules ? "${params.modules}".split(',') : []

/*************************************************
* IMPORT MODULES
**************************************************/
include { CREATE_DECOYS_FILE } from '../modules/salmon/decoys'
include { SALMON_INDEX       } from '../modules/salmon/index'
include { TX2GENE            } from '../modules/genomicfeatures/tx2gene'

/*************************************************
* WORKFLOW
**************************************************/
workflow CHECK_AND_PREPARE_INPUT {
    main:
        
        //---------------------------------------------
        // INITIALISE OUTPUT
        //---------------------------------------------
        ch_reads         = Channel.empty()
        ch_salmon_index  = Channel.empty()
        ch_genome        = Channel.empty()
        ch_transcriptome = Channel.empty()
        ch_gtf           = Channel.empty()
        ch_genemap       = Channel.empty()
        ch_versions      = Channel.empty()

        //---------------------------------------------
        // CHECK PARAMS AND DEFINE VARIABLES
        //---------------------------------------------
        def prefix  = ( params.prefix && (params.prefix  instanceof String)) ? params.prefix  : ''
        def suffix1 = (params.suffix1 && (params.suffix1 instanceof String)) ? params.suffix1 : ''
        def suffix2 = (params.suffix2 && (params.suffix2 instanceof String)) ? params.suffix2 : ''
        def refdir  = ( params.refdir && (params.refdir  instanceof String)) ? params.refdir  : "${params.outdir}/ref/${params.species}"

        def download_transcriptome = false
        def download_genome = false

        //---------------------------------------------
        // CHECK READS
        //---------------------------------------------
        //Check mandatory input is provided
        if (! params.input) { 
            exit 1, 'ERROR: `input` not specified!' 
        } else {
            file(params.input, checkIfExists: true)
        }

        //Check file extension
        if( !(supported_filext.contains(params.filext)) ){
            exit 1, "ERROR: file extension `${params.filext}` not supported"
        }

        //Create input channel
        ch_reads = create_channel_from_dir(params.input, params.filext, suffix1, suffix2, prefix)
 
        //Modules to run
        if( modules_to_run.size() > 0 ){
            if (! supported_modules.containsAll(modules_to_run) ){
                error "ERROR: provided pipeline modules not matching supported types. Available options: 'fastqc', 'quant', 'multiqc'. Example: --modules fastqc,quant"
            }
        } else {
            //error "ERROR: no pipeline module was selected to run. Available options: 'fastqc', 'quant', 'multiqc'. Example: --modules fastqc,quant"
            log.info "\nWARNING: no pipeline module was selected to run\n"
        }


        //---------------------------------------------
        // CHECK REFERENCE
        //---------------------------------------------
        //TRANSCRIPTOME
        if(params.transcriptome && file(params.transcriptome).exists()){
            log.info "User-provided transcriptome file found"
            ch_transcriptome = Channel.fromPath(params.transcriptome, type: 'file')
        } else {
            if(raw_reads_filext.contains(params.filext)){
                ch_transcriptome = Channel.empty()
                download_transcriptome = true
            } else {
                exit 1, "ERROR: aligned reads provided in input (file extension `${params.filext}`) but no reference transcriptome provided"
            }
        }
            
        //GENOME
        if(params.genome && file(params.genome).exists()){
            log.info "User-provided genome file found"
            ch_genome = Channel.fromPath(params.genome, type: 'file')
        } else {
            ch_genome = Channel.empty()
            download_genome = true
        }

        if(modules_to_run.contains('quant')){
            //SALMON INDEX
            if (! params.salmon_index) { 
                exit 1, 'ERROR: `salmon_index` not specified!' 
            } else {
                def salmon_index = file(params.salmon_index)//returns file system object
                if ( salmon_index.exists() && salmon_index.getExtension()=='gz') {
                    log.info "Compressed salmon index found: decompressing"
                    // Uncompress
                    GUNZIP(
                        salmon_index
                    )
                    // Update
                    salmon_index = GUNZIP.out.decompressed
                } 
                
                if ( salmon_index.exists() && salmon_index.isDirectory() && salmon_index.listFiles().size() > 0) {
                    log.info "Salmon index directory found and not empty"
                    ch_salmon_index = salmon_index
                } else {
                    log.info "Salmon index directory not found or empty: a new transcriptome index will be created"
                    
                    // Check transcriptome
                    if( download_transcriptome ){
                        //Download
                        log.info "User-provided transcriptome file not found: try to download it"
                        // Check species
                        if( !params.species ){
                            exit 1, "ERROR: `species` not specified!"
                        }
                        // Get URL
                        def url = getTranscriptomeUrl(params.species)
                        // Download file
                        transcriptomeFilePath = downloadFile(
                            url, 
                            refdir, 
                            null,
                            true
                        )
                        // Update channel
                        ch_transcriptome = Channel
                            .fromPath(transcriptomeFilePath, checkIfExists: true, type: 'file')
                    } 
                    
                    
                    // Check genome
                    if( download_genome ){
                        //Download
                        log.info "User-provided genome file not found: try to download it"
                        // Check species
                        if( !params.species ){
                            exit 1, "ERROR: `species` not specified!"
                        }
                        // Get URL
                        def url = getGenomeUrl(params.species)
                        // Download file
                        genomeFilePath = downloadFile(
                            url, 
                            refdir, 
                            null,
                            true
                        )
                        // Update channel
                        ch_genome = Channel
                            .fromPath(genomeFilePath, checkIfExists: true, type: 'file')
                    }

                    // Check decoys
                    if(params.decoys){
                        if((params.decoys instanceof String) && file(params.decoys).exists()){
                            log.info "Decoys file found"
                            ch_decoys = Channel.fromPath(params.decoys, type: 'file')
                        } else {
				            if(params.decoys instanceof String){
				                fn_decoys = file(params.decoys).name
				            } else {
					            fn_decoys = "decoys.txt"
				            }
                            log.info "Decoys file not found. A decoys file will be generated."
                            // Create decoys file.
			                CREATE_DECOYS_FILE(
                                ch_genome,
                                fn_decoys
                                // file(params.decoys).parent
                            )
                            ch_decoys = CREATE_DECOYS_FILE.out.decoys
                            // Update versions
                            ch_versions = ch_versions.mix(CREATE_DECOYS_FILE.out.versions) 
                        }
                    } else {
                        log.info "Decoys option is not selected. No decoys will be used in generating the index. This is not recommended."
                        ch_decoys = Channel.empty()
			        }

                    // Create index
                    SALMON_INDEX(
                        ch_transcriptome,
                        ch_genome,
                        ch_decoys.collect().ifEmpty([])
                    )

                    // Output channel
                    ch_salmon_index = SALMON_INDEX.out.index
                    
                    // Update versions
                    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)     
                }
            }

            //GTF
            if(params.gtf){
                if((params.gtf instanceof String) && file(params.gtf).exists()){
                    log.info "User-provided gtf file found"
                    ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true, type: 'file')
                } else {
                    //Download
                    log.info "User-provided gtf file not found: try to download it"
                    // Check species
                    if( !params.species ){
                        exit 1, "ERROR: `species` not specified!"
                    }
                    // Get URL
                    def url = getGtfUrl(params.species)
                    // Download file
                    gtfFilePath = downloadFile(
                        url, 
                        refdir, 
                        null,
                        true
                    )
                    // Update channel
                    ch_gtf = Channel
                        .fromPath(gtfFilePath, checkIfExists: true, type: 'file')
                }
            } else {
                ch_gtf = Channel.empty()
            }

            // GENEMAP
            if(params.genemap){
                if((params.genemap instanceof String) && file(params.genemap).exists()){
                    log.info "User-provided genemap file found"
                    ch_genemap = Channel.fromPath(params.genemap, checkIfExists: true, type: 'file')
                } else if(ch_gtf){
                    // Create from GTF
                    log.info "User-provided genemap file not found: try to create it from GTF file"
                    // Process GTF
                    TX2GENE(ch_gtf)
                    // Update channel
                    ch_genemap = TX2GENE.out.genemap
                    // Update versions
                    ch_versions = ch_versions.mix(TX2GENE.out.versions)
                } else {
                    ch_genemap = Channel.empty()
                }
            } else {
                ch_genemap = Channel.empty()
            }
        }




    emit:
        reads         = ch_reads            // channel: [ val(meta), [ reads1 ], [ reads2 ] ]
        salmon_index  = ch_salmon_index     // channel: /path/to/salmon/index
        genome        = ch_genome           // channel: /path/to/genome.fa.gz
        transcriptome = ch_transcriptome    // channel: /path/to/transcriptome.fa.gz
        gtf           = ch_gtf              // channel: /path/to/gtf.gtf
        genemap       = ch_genemap          // channel: /path/to/tx2gene.tsv
        versions      = ch_versions         // channel: versions.yml
}
