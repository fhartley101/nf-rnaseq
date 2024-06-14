/*  A Nextflow process block
*/
process CUTADAPT {
     /*********** DIRECTIVES ***********/
    //Set to true to forward process stdout to the current top stdout
    debug false
    //Associate process execution with custom label
    tag "Cutadapt on $meta.id"
    //Annotate process with identifier
    label "process_medium"
    //Process dependencies
    conda "bioconda::cutadapt=4.8-0"
    //To execute the script in a Singularity or Docker container
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.8--py39hf95cd2a_0' : 
        'quay.io/biocontainers/cutadapt:4.8--py39hf95cd2a_0'}"
    //Specifies where to publish output
    publishDir path: "${params.outdir}/cutadapt" , mode: "${params.publishdir_mode}", overwrite: true, followLinks: true
     
    /*********** INPUT ***********/
    input:
        tuple val(meta), path(reads1), path(reads2)
	    
    /*********** OUTPUT ***********/
    //Define output channels and assign identifiers
    output:
        tuple val(meta), path("${meta.id}_1.trimmed.${params.filext}"), path("${meta.id}_2.trimmed.${params.filext}", arity: '0..2'), emit: reads // Ideally this would be arity: '0..1' but this is not yet supported. When the nullable option is added to Nextflow, that would be the optimum output option for read 2
        path "versions.yml"                 , emit: versions
        path "${meta.id}_info.txt"          , emit: info
    
    /*********** SCRIPT ***********/
    script:
        def args = task.ext.args ?: ''

        // Determine input
        def files1 = reads1.join(' ')
        def files2 = reads2.join(' ')
	
	if(meta.single_end){
	
	log.info "Cutadapt running on single end data with: ${args}"
	
	"""
	#!/usr/bin/env bash
	
    cutadapt \\
        -o ${meta.id}_1.trimmed.${params.filext} \\
        ${args} \\
        ${files1} 1>> ${meta.id}_info.txt
	
	sed -i "2s/\$/ # ${meta.id}.fastq.gz/" ${meta.id}_info.txt
	
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    cutadapt: \$(echo \$(cutadapt --version))
	END_VERSIONS
        
	"""
	
	} else {
	
	log.info "Cutadapt running on paired end data with: ${args}"
	
	"""
	#!/usr/bin/env bash
	
    cutadapt \\
        -o ${meta.id}_1.trimmed.${params.filext} \\
	    -p ${meta.id}_2.trimmed.${params.filext} \\
        ${args} \\
	    ${files1} ${files2} 1>> ${meta.id}_info.txt

	sed -i "2s/\$/ # ${meta.id}.fastq.gz/" ${meta.id}_info.txt
	
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    cutadapt: \$(echo \$(cutadapt --version))
	END_VERSIONS
        """
	}
	 
}
