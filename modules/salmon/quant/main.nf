/*  A Nextflow process block
*/
process SALMON_QUANT {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "Salmon on $meta.id"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    conda "bioconda::salmon=1.10.0"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.0--h7e5ed60_0' : 
        'quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/salmon/quant" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(reads1), path(reads2)
        path index
        path transcriptome
        path genemap
        val libtype
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        tuple val(meta), path("${meta.id}") , emit: quant
        path "versions.yml"                 , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        
        // Check if the threads option is passed from config, use task.cpus otherwise
        def threads = (args.indexOf('-p ') > 1) ? '' : "-p ${task.cpus}"
        
        // Determine which mode to run based on the reads file extension (run mapping-based vs alignment-based)
        String[] mapping_mode_filext = ['fa', 'fasta', 'fa.gz', 'fasta.gz', 'fq', 'fastq',  'fastq.gz', 'fq.gz']
        String[] alignment_mode_filext = ['bam', 'sam']
        def salmon_mode;
        def reads = reads1.toList().toArray()
        def has_mapping_mode_ext;
        def has_alignment_mode_ext;
        // Check all files have extensions belonging to the same mode
        if(reads.size()>1){
            has_mapping_mode_ext   = Arrays.stream(reads).allMatch(ifilepath -> Arrays.stream(mapping_mode_filext).anyMatch(entry -> ifilepath.name.endsWith(entry)))
            has_alignment_mode_ext = Arrays.stream(reads).allMatch(ifilepath -> Arrays.stream(alignment_mode_filext).anyMatch(entry -> ifilepath.name.endsWith(entry)))
        } else {
            ifilepath = reads[0]
            has_mapping_mode_ext   = Arrays.stream(mapping_mode_filext).anyMatch(entry -> ifilepath.name.endsWith(entry))
            has_alignment_mode_ext = Arrays.stream(alignment_mode_filext).anyMatch(entry -> ifilepath.name.endsWith(entry))
        }
        
        if ( has_mapping_mode_ext ){
            if (index){
                salmon_mode = 'mapping-based'
            } else {
                error "Invalid input: 'index' must be provided for mapping-based mode."
            }
        } else if ( has_alignment_mode_ext ){
            if (transcriptome){
                salmon_mode = 'alignment-based'
            } else {
                error "Invalid input: 'transcriptome' must be provided for alignment-based mode."
            }
        } else {
            error "Invalid input: 'reads' must have a supported file extension (${mapping_mode_filext.join(', ')}, ${alignment_mode_filext.join(', ')})."
        }

        // Set input reads and reference
        def is_mapping_mode = (salmon_mode == 'mapping-based') ? true : false
        // log.info "Mapping mode: ${is_mapping_mode}"
        // Determine input
        def input_reads;
        def reference;
        def files1 = reads1.join(' ')
        def files2 = reads2.join(' ')
        if (is_mapping_mode){
            // Check if single-end or paired-end reads    
            input_reads = meta.single_end ? "--unmatedReads ${files1}" : "--mates1 ${files1} --mates2 ${files2}"
            // Set index
            reference = "--index ${index}"
        } else {
            input_reads = "-a ${files1}"
            reference = "-t ${transcriptome}"
        }
        // log.info "Reference: ${reference}"
        
        //Supported library types
        String[] supported_libtypes = [
            'I'  , 'O'  , 'M'  ,//inward,outward,matching (only for paired-end reads)
            'U'  ,              //unstranded
            'IU' , 'OU' , 'MU' ,//unstranded
            'IS' , 'OS' , 'MS' ,//stranded
            'ISF', 'OSF', 'MSF',//read 1 (or single-end read) comes from the forward strand
            'ISR', 'OSR', 'MSR',//read 1 (or single-end read) comes from the reverse strand
            'SF' , 'SR' ,       //stranded
            'A'                 //automatic
        ]
        //Check provided library type
        String libType = 'A'
        if (libtype && supported_libtypes.contains(libtype)){
            libType = libtype
        } else {
            log.info "Invalid library type specified (${libtype}), defaulting to auto-detection with '--libType=A'."
        }

        //Check genemap
        def geneMap = (genemap.name != 'NO_FILE') ? "--geneMap ${genemap}" : ''
        // log.info "Genemap: ${geneMap}"
        
        """
        #!/usr/bin/env bash
        
        salmon quant --libType=${libType} ${reference} \\
            ${input_reads} \\
            --output ${meta.id} \\
            ${threads} \\
            ${geneMap} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
}