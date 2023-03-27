process MULTIQC {
    /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Annotate process with identifier
    label "process_single"
    //Process dependencies
    conda "bioconda::multiqc=1.14"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' : 
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/multiqc" , mode: "${params.publishdir_mode}", overwrite: true
    
    /*********** INPUT ***********/
    input:
        path multiqc_config
        path software_versions
        path ("fastqc/*")
        path ("salmon/*")
    
    /*********** OUTPUT ***********/
    output:
        path "*multiqc_report.html", emit: report
        path "*_data"              , emit: data
        path "*_plots"             , optional:true, emit: plots
        path "versions.yml"        , emit: versions
        
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        def custom_config = multiqc_config ? "--config $multiqc_config" : ''

        """
        #!/usr/bin/env bash
        echo ${multiqc_config}
        multiqc \\
            -f \\
            ${args} \\
            ${multiqc_config} \\
            .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
}