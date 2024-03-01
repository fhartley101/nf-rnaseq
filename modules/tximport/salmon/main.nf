/*  A Nextflow process block
*/
process TXIMPORT_SALMON {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "TXIMPORT_SALMON"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    conda "bioconda::bioconductor-tximeta=1.12.0"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.16.0--r42hdfd78af_0' : 
        'quay.io/biocontainers/bioconductor-tximeta:1.16.0--r42hdfd78af_0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/salmon/tximport" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path ("salmon/*")
        path  tx2gene
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path "*gene_tpm.tsv"                 , emit: tpm_gene
        path "*gene_counts.tsv"              , emit: counts_gene
        path "*gene_length_scaled_counts.tsv", emit: counts_gene_length_scaled
        path "*gene_scaled_counts.tsv"       , emit: counts_gene_scaled
        path "*transcript_tpm.tsv"           , emit: tpm_transcript
        path "*transcript_counts.tsv"        , emit: counts_transcript
        path "*tximport_gene_summary.rds"    , emit: tximport_gene_summary
        path "versions.yml"                  , emit: versions
    
    /*********** SCRIPT ***********/
    script: // This script is bundled with the pipeline, in /bin/

        """
        #!/usr/bin/env bash
        tximport_salmon.R salmon ${tx2gene}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            bioconductor-tximport: \$(Rscript -e "library(tximport); cat(as.character(packageVersion('tximport')))")
        END_VERSIONS
        """
}
