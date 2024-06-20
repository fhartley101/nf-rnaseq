/*  A Nextflow process block
*/
process CREATE_DECOYS_FILE {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug true
    //Associate process execution with custom label
    tag "decoys"
    //Annotate process with identifier
    label "process_low"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/salmon/" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path genome
        val decoys
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path decoys         , emit: decoys 
        path "versions.yml" , emit: versions
        
    /*********** SCRIPT ***********/
    script:
        // Output 
        def decoys = decoys ?: "decoys.txt"

        // Log info
        log.info "decoys filename ${decoys}"
        log.info "decoys outdir ${params.refdir}"

        // File extension 
        def dfilext = "${genome.extension}"
        // log.info "extension ${dfilext}"

        def is_compressed;
        if ( genome.name.endsWith('.gz') ){
            log.info "Genome is compressed gz file"
            is_compressed=true
        } else if ( genome.name.endsWith('.fa') || genome.name.endsWith('.fasta')){
            log.info "Genome is uncompressed FASTA file"
            is_compressed=false
        } else {
            error "ERROR: genome file extension (${dfilext}) not supported"
        }

        log.info "Creating a new decoys file" 

        if( is_compressed ){
            """
            #!/usr/bin/env bash
            
            echo "Decompressing the reference genome file"
            ugfname=`basename "${genome}" .gz`      
            gunzip -d -k -c "${genome}" > "${genome.baseName}"
            echo "Extracting the names of the genome targets"
            grep "^>" < "\${ugfname}" | cut -d " " -f 1 | cut -f 1 > "${decoys}"
            
            echo "Editing the decoys file"
            sed -i.bak -e 's/>//g' "${decoys}"

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                gunzip: \$(gunzip --version 2>&1 | head -n 1 | sed -e 's/^.*gzip) //g' )
                cut: \$(cut --version 2>&1 | head -n 1 | sed -e 's/^.*coreutils) //g' )
                grep: \$(grep --version 2>&1 | head -n 1 | sed -e 's/^.*grep) //g' ) 
            END_VERSIONS

            """
        } else {
            """
            #!/usr/bin/env bash
            
            echo "Extracting the names of the genome targets"
            grep "^>" < "${genome}" | cut -d " " -f 1 | cut -f 1 > "${decoys}"

            echo "Editing the decoys file"
            sed -i.bak -e 's/>//g' "${decoys}"

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cut: \$(cut --version 2>&1 | head -n 1 | sed -e 's/^.*coreutils) //g')
                grep: \$(grep --version 2>&1 | head -n 1 | sed -e 's/^.*grep) //g' )
            END_VERSIONS
            """
        }
}

