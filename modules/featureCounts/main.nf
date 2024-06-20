/*  A Nextflow process block
*/
process FEATURE_COUNTS {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "featureCounts on $meta.id"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    conda "bioconda::subread=2.0.6-0"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_2' :
        'quay.io/biocontainers/subread:2.0.6--he4a0461_2'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/featureCounts", mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(bam)
        path gtf
        val stranded
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        tuple val(meta), path ("${meta.id}.featurecounts") , emit: featurecounts
        path "${meta.id}.featurecounts"                    , emit: counts
        path "${meta.id}.featurecounts.summary"            , emit: summary
        path "versions.yml"                                , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        def threads = (args.indexOf('-T ') > 1) ? '' : "-T ${task.cpus}"

        def paired_end = (!meta.single_end) ? "-p" : ""

        def gtffilename = "${gtf.name}"


        """
        #!/usr/bin/env bash

        featureCounts ${threads} \
            ${args} \
            -a ${gtffilename} \
            -o ${meta.id}.featurecounts \
            ${paired_end} \
            -s ${stranded} \
            ${bam}

        cat <<-END_VERSIONS > versions.yml
        "$task.process":
            featureCounts: \$(echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" )
        END_VERSIONS

        """
}