process SOFTWARE_VERSIONS {
    /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Annotate process with identifier
    label "process_single"
    //Process dependencies
    conda "conda-forge::r-base=4.2.1"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'quay.io/biocontainers/r-base:4.2.1'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/pipeline_info" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path versions

    /*********** OUTPUT ***********/
    output:   
        path "software_versions.yml"        , emit: yml
        path "software_versions_mqc.yml"    , emit: multiqc_yml
        path "software_versions.html"       , emit: html

    /*********** SCRIPT ***********/

    script: 
        def args = task.ext.args ?: ''
        """
        #!/usr/bin/env bash
        software_versions.R ${versions}
        """
}