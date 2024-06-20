/*  A Nextflow process block
*/
process SALMON_INDEX {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "salmon_index"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    conda "bioconda::salmon=1.10.0"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.0--h7e5ed60_0' : 
        'quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/salmon", mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path transcriptome
        path genome
        path decoys
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path "./index"     , emit: index
        path "versions.yml" , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        def threads = (args.indexOf('-p ') > 1) ? '' : "-p ${task.cpus}"
        def tfilename = "${transcriptome.name}"
        def gfilename = "${genome.name}"
        def gentrome;
        // Log
        // log.info "Genome ${gfilename}"
        // log.info "Transcriptome ${tfilename}"
        
        // Check extension
        if( (tfilename.endsWith('.fa.gz') || tfilename.endsWith('.fasta.gz')) &&  
            (gfilename.endsWith('.fa.gz') || gfilename.endsWith('.fasta.gz'))){
            gentrome = 'gentrome.fa.gz'
            log.info "Genome and transcriptome are both compressed FASTA files"
        } else if( (tfilename.endsWith('.fa') || tfilename.endsWith('.fasta')) &&
            (gfilename.endsWith('.fa') || gfilename.endsWith('.fasta'))){
            gentrome = 'gentrome.fa'
            log.info "Genome and transcriptome are both uncompressed FASTA files"
        } else {
            error "ERROR: genome and transcriptome extensions do not match"
        }

        // Check decoys
        def idecoys = decoys ? "--d ${decoys}" : ''

        """
        #!/usr/bin/env bash

        cat ${transcriptome} ${genome} > ${gentrome}

        salmon index \
            ${threads} \
            -t ${gentrome} \
            ${idecoys} \
            ${args} \
            -i index

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
}