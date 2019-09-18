#!/usr/bin/env Rscript

cat(R.version$version.string, "\n")

# load arguments ---------------------------------------------------------
args <- commandArgs(TRUE)
inp.abundances.path  <- args[[1]] # OTU table file 
inp.metadata.path <- args[[2]] # Metadata file 
diffsource <- args[[3]] # diff source (0,1)
output <- args[[4]] # to write to temp file

# load data ---------------------------------------------------------------

# table import 
if(!file.exists(inp.abundances.path)) {
  errQuit("Input count table file does not exist.")
} else {
  table <- Load_CountMatrix(CountMatrix_path = inp.abundances.path)
}
# map import
if(!file.exists(inp.metadata.path)) {
  errQuit("Input mapping file does not exist.")
} else {
  map <- Load_metadata(metadata_path = inp.metadata.path)
}

# load libraries ----------------------------------------------------------
suppressWarnings(library(FEAST))

# analysis ----------------------------------------------------------------
proportions <- FEAST(C = table,
                     metadata = meta,
                     different_sources_flag = diffsource)

# output ----------------------------------------------------------------
proportions <- as.data.frame(proportions)
write.csv(proportions, file=output)