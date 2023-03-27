/*  A Nextflow process block
*/
process GUNZIP {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug true
    //Annotate process with identifier
    label "process_single"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/decompressed" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path compressed
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path decompressed   , emit: decompressed 
        path "versions.yml" , emit: versions
        
    /*********** SCRIPT ***********/
    script:
        // Output
        def decompressed = compressed.baseName 

        if ( !compressed.name.endsWith('.gz') ){
            error 1, "Input file is not a compressed `gz` file"
        } 

        """
        #!/usr/bin/env bash

        gunzip -d -k -c "${compressed}" > "${decompressed}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gunzip: \$(gunzip --version 2>&1) 
        END_VERSIONS

        """
}

