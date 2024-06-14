/*  A Nextflow process block
*/
process STAR_INDEX {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "star_index"
    //Annotate process with identifier
    label "process_high"
    //Process dependencies
    conda "bioconda::star=2.7.11b-1"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/star:2.7.8a--0' :
        'quay.io/biocontainers/star:2.7.8a--0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/star", mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
    
    /*********** INPUT ***********/
    input:
        path gtf
        path genome
    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        path "index/"        , emit: index
        path "versions.yml"       , emit: versions
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''
        def threads = (args.indexOf('--runThreadN ') > 1) ? '' : "--runThreadN ${task.cpus}"

        def gtffilename = "${gtf.name}"
        def gfilename = "${genome.name}"

        // Remove the gz extension if the file is gzipped. This does not affect uncompressed file names
        def gtfbasename = gtffilename.replace(".gz", "")
        def gbasename = gfilename.replace(".gz", "")

        """
        #!/usr/bin/env bash
       
        mkdir index

        gunzip ${gtffilename}
        gunzip ${gfilename}

        STAR --runMode genomeGenerate \
            ${threads} \
            --sjdbGTFfile ${gtfbasename} \
            --genomeFastaFiles ${gbasename} \
            ${args} \
            --genomeDir index/

        rm ${gtfbasename}
        rm ${gbasename}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(echo \$(STAR --version))
        END_VERSIONS

        """

}