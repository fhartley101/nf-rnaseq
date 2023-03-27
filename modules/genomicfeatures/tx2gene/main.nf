process TX2GENE {
    /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "TX2GENE"
    //Annotate process with identifier
    label "process_single"
    //Process dependencies
    conda "bioconda::bioconductor-genomicfeatures=1.50.2"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bioconductor-genomicfeatures:1.50.2--r42hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-genomicfeatures:1.50.2--r42hdfd78af_0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/ref/${params.species}/" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path gtf

    /*********** OUTPUT ***********/
    output:   
        path 'tx2gene.tsv'  , emit: genemap
        path "versions.yml" , emit: versions

    /*********** SCRIPT ***********/
    script: // This script is bundled with the pipeline, in /bin/
        """
        #!/usr/bin/env bash
        tx2gene.R ${gtf}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            bioconductor-genomicranges: \$(Rscript -e "library(GenomicFeatures); cat(as.character(packageVersion('GenomicFeatures')))")
        END_VERSIONS
        """
}