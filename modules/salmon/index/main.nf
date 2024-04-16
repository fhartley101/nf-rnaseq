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
    publishDir path: "${params.outdir}/index/${params.species}", mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path transcriptome
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path "./salmon"     , emit: index
        path "versions.yml" , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        def threads = (args.indexOf('-p ') > 1) ? '' : "-p ${task.cpus}"
        def tfilename = "${transcriptome.name}"
        // Log
        // log.info "Transcriptome ${tfilename}"
        
        // Check extension
        if (tfilename.endsWith('.fa.gz') || tfilename.endsWith('.fasta.gz')){
            log.info "Transcriptome is a compressed FASTA file"
        } else (tfilename.endsWith('.fa') || tfilename.endsWith('.fasta')){
            log.info "Transcriptome is an uncompressed FASTA file"
        }

        """
        #!/usr/bin/env bash

        salmon index \
            ${threads} \
            -t ${transcriptome} \
            ${args} \
            -i salmon

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
}
