/*  A Nextflow process block
*/
process GENECOUNTS {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "GENECOUNTS"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    //conda "bioconda::bioconductor-tximeta=1.12.0"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' : 
        'quay.io/biocontainers/pandas:1.5.2'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/featureCounts" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path ("featureCounts/*")
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path "genecounts.csv"           , emit: genecounts
        path "versions.yml"             , emit: versions
    
    /*********** SCRIPT ***********/
    script: // This script is bundled with the pipeline, in /bin/
        
        """
        #!/usr/bin/env bash
        genecounts.py featureCounts

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python3 --version) | sed -e "s/Python //g") 
            pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
        END_VERSIONS

        """
}
