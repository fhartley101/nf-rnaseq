process CAT_SINGLE_READS {
    /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "$meta.id"
    //Annotate process with identifier
    label "process_single"
    
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/input/merged" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(reads1), path(reads2)
        val verbose

    /*********** OUTPUT ***********/
    output:   
        tuple val(meta), path("${meta.id}_1.*" ), path(''), emit: reads
        path "versions.yml"                               , emit: versions

    /*********** SCRIPT ***********/
    script:
        // Extract file name and extension 
        // NB Input is supposed to have same file extension, so we can use first element
        def filename = "${reads1[0].name}"
        def filext = "${reads1[0].extension}"

        // Determine if is compressed file
        def is_compressed = filename.endsWith('.gz')
        if(is_compressed){
            def bName = "${reads1[0].baseName}"
            filext = bName.substring(bName.lastIndexOf(".") + 1)
            filext = "${filext}.gz"
        }
        if(verbose) { log.info "File ext: ${filext}" }
        // Convert to string
        def reads1List = (reads1 instanceof List || reads1 instanceof Arrays) ? reads1.collect{ it.toString() } : [reads1.toString()]

        if(verbose) { log.info "readList: ${reads1List}" }
        if (reads1List.size >= 1) {
            if (meta.single_end) {
                if(is_compressed){
                    """
                    #!/usr/bin/env bash
                    gunzip -c ${reads1List.join(' ')} | gzip > ${meta.id}_1.${filext}

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        gunzip: \$(gunzip --version 2>&1) 
                    END_VERSIONS
                    """
                } else {
                    """
                    #!/usr/bin/env bash
                    cat ${reads1List.join(' ')} > ${meta.id}_1.${filext}

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        cat: echo \$(cat --version 2>&1) | sed 's/^.*(cat) //; s/ Copyright.*\$//'
                    END_VERSIONS
                    """ 
                }

            } 
        } 
}

process CAT_PAIRED_READS {
    /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "$meta.id"
    //Annotate process with identifier
    label "process_single"
    
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/input/merged" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(reads1), path(reads2)
        val verbose

    /*********** OUTPUT ***********/
    output:   
        tuple val(meta), path("${meta.id}_1.*" ), path("${meta.id}_2.*" ), emit: reads
        path "versions.yml"                                              , emit: versions

    /*********** SCRIPT ***********/
    script:
        // Extract file name and extension 
        // NB Input is supposed to have same file extension, so we can use first element
        def filename = "${reads1[0].name}"
        def filext = "${reads1[0].extension}"

        // Determine if is compressed file 
        def is_compressed = filename.endsWith('.gz')
        if(is_compressed){
            def bName = "${reads1[0].baseName}"
            filext = bName.substring(bName.lastIndexOf(".") + 1)
            filext = "${filext}.gz"
        }
        if(verbose) { log.info "File ext: ${filext}" }
        // Convert to string
        def reads1List = (reads1 instanceof List || reads1 instanceof Arrays) ? reads1.collect{ it.toString() } : [reads1.toString()]
        def reads2List = (reads2 instanceof List || reads2 instanceof Arrays) ? reads2.collect{ it.toString() } : [reads2.toString()]

        if(verbose) { log.info "readList: ${reads1List}" }
        if (reads1List.size >= 1) {
            if (meta.single_end) {
                if(is_compressed){
                    """
                    #!/usr/bin/env bash
                    gunzip -c ${reads1List.join(' ')} | gzip > ${meta.id}_1.${filext}

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        gunzip: \$(gunzip --version 2>&1) 
                    END_VERSIONS
                    """
                } else {
                    """
                    #!/usr/bin/env bash
                    cat ${reads1List.join(' ')} > ${meta.id}_1.${filext}

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        cat: echo \$(cat --version 2>&1) | sed 's/^.*(cat) //; s/ Copyright.*\$//'
                    END_VERSIONS
                    """ 
                }

            } else {
                if(is_compressed){
                    """
                    #!/usr/bin/env bash
                    gunzip -c ${reads1List.join(' ')} | gzip > ${meta.id}_1.${filext}
                    gunzip -c ${reads2List.join(' ')} | gzip > ${meta.id}_2.${filext}

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        gunzip: \$(gunzip --version 2>&1) 
                    END_VERSIONS
                    """
                } else {
                    """
                    #!/usr/bin/env bash
                    cat ${reads1List.join(' ')} > ${meta.id}_1.${filext}
                    cat ${reads2List.join(' ')} > ${meta.id}_2.${filext}

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        cat: echo \$(cat --version 2>&1) | sed 's/^.*(cat) //; s/ Copyright.*\$//'
                    END_VERSIONS
                    """ 
                }
            }
        } 
}