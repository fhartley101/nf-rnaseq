#!/usr/bin/env Rscript

# Command Args -----------------------------------------------------------------
# Extract command line arguments
args = commandArgs(trailingOnly = TRUE)
# Check usage
if (length(args) < 1) {
  stop("Usage: software_versions.r <versions.yml>", call. = FALSE)
}
# Path to YAML file containing the versions of used software
versions = args[1]

# Read -------------------------------------------------------------------------
# File path
fp = file.path(versions)

# Parse file
data = readLines(
  con = fp
)

# Check if any illegal row -----------------------------------------------------
# Split each row
datasplit = strsplit(x = data, split = ":")
# Select rows with >2 elements
rows_to_check = which(lapply(X = datasplit, FUN = length)>2)
# Check row
row_to_check = rows_to_check[1]
for (row_to_check in rows_to_check){
  # Get original row
  tmp = data[[row_to_check]]
  # Check if string contains illegal/unrecognised option
  is_bad = any(grepl(pattern = 'illegal option|unrecognized option|unrecognised option', x = tmp))
  # If bad, change string
  if(isTRUE(is_bad)){
    data[[row_to_check]] = paste0(datasplit[[row_to_check]][1], ": NA")
  }
}

# Map --------------------------------------------------------------------------
# List
datal = as.list(data)
outl  = list()

# Fill list
mapi = 0
for (i in seq_len(length(datal))) {
  el = datal[[i]]
  # Check is not empty
  if(isTRUE(el=='')){next;}
  # Extract first char
  chr1 = substr(x = el, start = 1, stop = 2)
  # Check if is key or element
  is_el = any(grepl(pattern = "[ \t]", x = el))
  if (isFALSE(is_el)) {
    # Clean string
    module_name = gsub(pattern = "[ \t\r\n\"]|:$", replacement = "", x = el)
    # Is module already present?
    has_module = module_name %in% names(outl)
    # Check
    if (isFALSE(has_module)) {
      # Add new module
      outl[[module_name]] = list()
      # Update index
      mapi = mapi + 1
    } else {
      # Message
      logMsg = paste0("Same module `",module_name,"` reported more software version files. It might be caused by an issue in the pipeline.")
      message(logMsg)
    }
    # Clean
    rm(module_name, has_module)
  } else {
    # Split
    trow = unlist(strsplit(x = el, split = ":"))
    # Clean strings
    trow = trimws(trow)
    # Get software info
    sw_name = trow[1]
    sw_version = trow[2]
    # Get module
    module = outl[[mapi]]
    # Has software name?
    has_sw_name = sw_name %in% names(module)
    # Check
    if (isFALSE(has_sw_name)) {
      # Add software version
      module = c(module, stats::setNames(object = list(sw_version), nm = sw_name))
      # Store
      outl[[mapi]] = module
    } else {
      # Get old sw version
      old_sw_version = module[[sw_name]]
      # Is same version?
      is_same_version = old_sw_version == sw_version
      # Check
      if (isFALSE(is_same_version)) {
        logMsg = paste0("Different software versions of `", sw_name, "` used across modules. Check the pipeline.")
        message(logMsg)
      }
    }
    # Clean
    rm(trow, sw_name, sw_version, module)
  }
}

# List to string as YAML -------------------------------------------------------

# Initialise
outyaml = ''
# Number of modules
nmodules = length(outl)
# Loop over modules
for(i in seq_len(nmodules)){
  # Module
  module = outl[[i]]
  # Module name
  module_name = names(outl)[i]
  # Write module name
  outyaml = paste0(outyaml, module_name, ":\n")
  # Loop over software
  for(software_name in names(module)){
    # Software version
    software_version = module[[software_name]]
    # Write software
    outyaml = paste0(
      outyaml,
      paste0("  ", software_name, ": ", software_version),
      "\n")
    # Clean
    rm(software_version)
  }
  # Clean
  rm(module, module_name)
}
# Clean
rm(i, nmodules)

# Data frame -------------------------------------------------------------------
# Create output dataframe
outdf = data.frame(
  module = character(),
  software = character(),
  version = character(),
  stringsAsFactors = F)

# Fill data frame
for (module_name in names(outl)){
  # Module
  module = outl[[module_name]]
  # Loop over software
  for (software_name in names(module)){
    # Version
    software_version = module[[software_name]]
    # Store
    outdf = rbind(
      outdf,
      data.frame(module_name, software_name, software_version, stringsAsFactors = F)
    )
  }
}

# HTML table -------------------------------------------------------------------
# Set HTML document
htmltable = "<!DOCTYPE html>\n"

# Create HTML table
htmltable = paste('<table class="table" style="width:100%" id="sw-versions">')

## Head ------------------------------------------------------------------------
htmltable = paste(htmltable, "<thead>")
# Fill head
for(cname in colnames(outdf)){
  htmltable = paste0(
    htmltable,
    "<th>",
    cname,
    "</th>"
  )
}
# End Head
htmltable = paste(htmltable, "</thead>")

## Body ------------------------------------------------------------------------
htmltable = paste(htmltable, "<tbody>")
for(irow in seq_len(nrow(outdf))){
  htmltable = paste(htmltable, "<tr>")
  for(icol in seq_len(ncol(outdf))){
    htmltable = paste0(
      htmltable,
      "<td>",
      outdf[irow, icol],
      "</td>"
    )
  }
  htmltable = paste(htmltable, "</tr>")
}

## End -------------------------------------------------------------------------
# End Body
htmltable = paste(htmltable, "</tbody>")
# End Table
htmltable = paste(htmltable, "</table>")
htmltable = paste(htmltable, "\n")

# MultiQC ----------------------------------------------------------------------
versions_multiqc = list(
  "id"          = "software_versions",
  "section_name"= "Software Versions",
  "description" = "Software versions collected at run time from the software output.",
  "plot_type"   = "html",
  "data"        = htmltable
)

# Initialise
multiqc_outyaml = ''
# Number of keys
nkeys = length(versions_multiqc)
# Loop over keys
for(i in seq_len(nkeys)){
  # Key name
  yaml_key = names(versions_multiqc)[i]
  # Key value
  yaml_value = versions_multiqc[[i]]
  # Write key:value
  multiqc_outyaml = paste0(
    multiqc_outyaml,
    yaml_key, ": ", yaml_value,
    "\n")
}
# Clean
rm(i, yaml_key, yaml_value)

# Save -------------------------------------------------------------------------
# Write YAML file
write.table(
  x         = outyaml,
  file      = "software_versions.yml",
  quote     = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

# Write HTML file
write.table(
  x         = htmltable,
  file      = "software_versions.html",
  quote     = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

# Write MultiQC YAML file
# MultiQC-specific data file filename must end in *_mqc
write.table(
  x         = multiqc_outyaml,
  file      = "software_versions_mqc.yml",
  quote     = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
