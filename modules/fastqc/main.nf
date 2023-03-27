process FASTQC {
    /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Annotate process with identifier
    label "process_low_medium"
    //Associate process execution with custom label
    tag "$meta.id"
    //Process dependencies
    conda "bioconda::fastqc=0.11.9"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1' : 
        'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/fastqc" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(reads1), path(reads2)
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        tuple val(meta), path("${meta.id}*/*.html"), emit: html
        tuple val(meta), path("${meta.id}*/*.zip") , emit: zip
        path  "versions.yml"                       , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        // println args
        // echo "\nFASTQC:   ${meta}\n"
        // echo "\nFASTQC:   ${reads}\n"
        def files1 = reads1.join(' ')
        def files2 = reads1.join(' ')

        if(meta.single_end){
            """
            #!/usr/bin/env bash
            
            mkdir ${meta.id}

            fastqc \\
                -o ${meta.id} \\
                --threads $task.cpus \\
                ${args} \\
                ${files1}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
            END_VERSIONS
            """
        } else {
            """
            #!/usr/bin/env bash
            
            mkdir ${meta.id}_1
            mkdir ${meta.id}_2

            fastqc \\
                -o ${meta.id}_1 \\
                --threads $task.cpus \\
                ${args} \\
                ${files1}
            
            fastqc \\
                -o ${meta.id}_2 \\
                --threads $task.cpus \\
                ${args} \\
                ${files2}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
            END_VERSIONS
            """
        }

}