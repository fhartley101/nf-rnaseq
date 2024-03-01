#!/usr/bin/env Rscript

library(tximport)

# Command Args -----------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: salmon_tximport.r <salmon_outdir> <tx2gene>", call.=FALSE)
}

# Root directory containing sub-directories, each being an output of `salmon quant`
salmon_quant_dirpath = args[1]
# tx2gene file path
tx2gene_filepath     = args[2]

# Setup ------------------------------------------------------------------------

# Salmon quant files
files = list.files(
  path         = salmon_quant_dirpath,
  pattern      = "quant.sf",
  all.files    = FALSE,
  full.names   = TRUE,
  recursive    = TRUE,
  ignore.case  = FALSE,
  include.dirs = FALSE
)

# Sample names
sample_names = basename(path = dirname(path = files))

# Transcript to gene map
tx2gene = read.table(
  file             = tx2gene_filepath,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)
colnames(tx2gene) = c("TXNAME", "GENEID")

# Check if transcript ids have '.' character
ignoreTxVersion = !grepl(pattern = '.', x = tx2gene$TXNAME[1], fixed = T)

# Import the data --------------------------------------------------------------
txi <- tximport::tximport(
  files   = files,
  type    = "salmon",
  tx2gene = tx2gene,
  txOut   = TRUE
)

colnames(txi$abundance) = sample_names
colnames(txi$counts)    = sample_names

# Summarise data ---------------------------------------------------------------
gene_s = tximport::summarizeToGene(
  object              = txi,
  tx2gene             = tx2gene,
  varReduce           = FALSE,
  ignoreTxVersion     = ignoreTxVersion,
  ignoreAfterBar      = FALSE,
  countsFromAbundance = "no"
)
gene_stpm = tximport::summarizeToGene(
  object              = txi,
  tx2gene             = tx2gene,
  varReduce           = FALSE,
  ignoreTxVersion     = ignoreTxVersion,
  ignoreAfterBar      = FALSE,
  countsFromAbundance = "scaledTPM"
)
gene_lstpm = tximport::summarizeToGene(
  object              = txi,
  tx2gene             = tx2gene,
  varReduce           = FALSE,
  ignoreTxVersion     = ignoreTxVersion,
  ignoreAfterBar      = FALSE,
  countsFromAbundance = "lengthScaledTPM"
)

# Save data --------------------------------------------------------------------

## Transcript-level ------------------------------------------------------------
write.table(
  x         = txi$abundance,
  file      = file.path("transcript_tpm.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  x         = txi$counts,
  file      = file.path("transcript_counts.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)

## Gene-level ------------------------------------------------------------------
# Summary at gene-level
## Save file for downstream differential gene expression analysis
saveRDS(object = gene_s, file = file.path("tximport_gene_summary.rds"))
## Save gene-level abundance and counts
write.table(
  x         = gene_s$abundance,
  file      = file.path("gene_tpm.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  x         = gene_s$counts,
  file      = file.path("gene_counts.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)
# Summary at gene-level scaled up to library size
write.table(
  x         = gene_stpm$abundance,
  file      = file.path("gene_scaled_tpm.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  x         = gene_stpm$counts,
  file      = file.path("gene_scaled_counts.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)
# Summary at gene-level scaled using the average transcript length over samples and then the library size
write.table(
  x         = gene_lstpm$abundance,
  file      = file.path("gene_length_scaled_tpm.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  x         = gene_lstpm$counts,
  file      = file.path("gene_length_scaled_counts.tsv"),
  quote     = FALSE,
  sep       = "\t",
  col.names = TRUE,
  row.names = TRUE
)

# Session info -----------------------------------------------------------------
# Print sessioninfo to standard out
citation("tximport")
sessionInfo()