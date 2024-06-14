/*  A Nextflow process block
*/
process STAR_QUANT {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "star_index"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    conda "bioconda::star=2.7.11b-1"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/star:2.7.8a--0' :
        'quay.io/biocontainers/star:2.7.8a--0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/star/quant", mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(reads1), path(reads2)
        path index
        path gtf
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path "${meta.id}_Aligned.out.bam"        , emit: bam
        path "${meta.id}_Log.final.out"          , emit: log
        path "versions.yml"                      , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        def threads = (args.indexOf('--runThreadN ') > 1) ? '' : "--runThreadN ${task.cpus}"

        def files1 = reads1.join(' ')
        def files2 = reads2.join(' ')
                        
        // Remove the gz extension if the file is gzipped. This does not affect uncompressed file names
        def gtffilename = "${gtf.name}"        
        def gtffilename_trim = gtffilename.replace(".gz", "")

        """
        #!/usr/bin/env bash

        gunzip ${gtffilename}

        STAR ${threads} \
            --genomeDir ${index} \
            --sjdbGTFfile ${gtffilename_trim} \
            --readFilesIn ${files1} ${files2} \
            --outFileNamePrefix ${meta.id}_ \
            --readFilesCommand zcat \
            ${args}

        rm Homo_sapiens.GRCh38.112.chr.gtf

        cat <<-END_VERSIONS > versions.yml
        $task.process:
            star: \$(echo \$(STAR --version))
        END_VERSIONS

        """
}