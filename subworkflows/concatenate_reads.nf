#!/usr/bin/env nextflow

/*************************************************
* DSL2
**************************************************/
//enable DSL2 syntax
nextflow.enable.dsl=2

/*************************************************
* IMPORT MODULES
**************************************************/
include { CAT_SINGLE_READS   } from '../modules/utils'
include { CAT_PAIRED_READS   } from '../modules/utils'

/*************************************************
* WORKFLOW
**************************************************/
workflow CONCATENATE_READS {
    take:
        reads           // channel: [ val(meta), [ path(reads1) ] [ path(reads2) ] ]
    
    main:
        
        //---------------------------------------------
        // INITIALISE OUTPUT
        //---------------------------------------------
        ch_reads         = Channel.empty()
        ch_versions      = Channel.empty()

        //---------------------------------------------
        // CONCATENATE READS
        //---------------------------------------------
        
        // Split input by number of elements (i.e., if multi-file)
        reads.branch {
            meta, reads1, reads2 ->
                single  : reads1.size() == 1
                    return [ meta, reads1.flatten(), reads2.flatten() ]
                multiple: reads1.size() > 1
                    return [ meta, reads1.flatten(), reads2.flatten() ]
        }.set { ch_reads }
        
        // ch_reads.single.view(it -> "[SINGLE-FILE] $it")
        // ch_reads.multiple.view(it -> "[MULTI-FILE] $it")

        // Run only if multi-file reads exist
        if( ch_reads.multiple ) {
            // Split into single and paired end reads
            ch_reads.multiple.branch {
                meta, reads1, reads2 ->
                    single  : meta.single_end == true 
                        return [ meta, reads1.flatten(), reads2.flatten() ]
                    paired: meta.single_end == false 
                        return [ meta, reads1.flatten(), reads2.flatten() ]
            }.set { ch_multiple_reads }

            // ch_multiple_reads.single.view(it -> "[SINGLE READS] $it")
            // ch_multiple_reads.paired.view(it -> "[PAIRED READS] $it")

            // Concatenate reads
            ch_merged_single_reads = CAT_SINGLE_READS (
                ch_multiple_reads.single,
                Channel.value(false)
            )
            ch_versions = ch_versions.mix(ch_merged_single_reads.versions.first())
            // ch_merged_single_reads.reads.view(it -> "[MERGED SINGLE READS] $it")
            ch_merged_paired_reads = CAT_PAIRED_READS (
                ch_multiple_reads.paired,
                Channel.value(false)
            )
            ch_versions = ch_versions.mix(ch_merged_paired_reads.versions.first())
            // ch_merged_paired_reads.reads.view(it -> "[MERGED PAIRED READS] $it")

            // Merged the concatenated reads
            ch_multiple_reads = ch_merged_single_reads
                .reads
                .mix(ch_merged_paired_reads.reads)
            // ch_multiple_reads.view(it -> "[MIXED CONCATENATED MULTI-FILE READS] $it")

            // Merge the concatenated reads with single-file reads
            ch_reads = ch_multiple_reads
                .mix(ch_reads.single)
            // ch_reads.view(it -> "[MIXED SINGLE AND MULTI READS] $it")

            // Make sure everything is in the correct form
            ch_reads = ch_reads.map {
                meta,reads1,reads2 ->
                reads1 = (reads1 instanceof List || reads1 instanceof Arrays) ? reads1 : [reads1]
                if(reads2){
                    reads2 = (reads2 instanceof List || reads2 instanceof Arrays) ? reads2 : [reads2]
                } else {
                    reads2 = []
                }
                
                return [meta, reads1, reads2]
            }
            ch_reads.view(it -> "[FINAL READS] $it")
        } else {
            // Return input data
            ch_reads = reads
        }
        
    
    emit:
        reads         = ch_reads            // channel: [ val(meta), [ reads1 ], [ reads2 ] ]
        versions      = ch_versions         // channel: versions.yml
}