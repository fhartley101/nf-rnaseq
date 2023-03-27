#!/usr/bin/env Rscript

library(GenomicFeatures)

# Command Args -----------------------------------------------------------------

# Extract command line arguments
args = commandArgs(trailingOnly=TRUE)

# Check usage
if (length(args) < 1) {
    stop("Usage: tx2gene.r <gtf>", call.=FALSE)
}

# Path to GTF file
gtf <- file.path(args[1])

# Transcript-To-Gene -----------------------------------------------------------

# Build a transcript database from GTF file
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)

# Get the names of all transcripts by using the keys function
k <- AnnotationDbi::keys(txdb, keytype="TXNAME")

# Return a dataframe
tx_map <- AnnotationDbi::select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")

# Save data --------------------------------------------------------------------

# Save
write.table(
    x         = tx_map, 
    file      = "tx2gene.tsv", 
    sep       = "\t", 
    quote     = FALSE, 
    row.names = FALSE
)

# Session info -----------------------------------------------------------------

# Print session info to standard out
citation("GenomicFeatures")
sessionInfo()