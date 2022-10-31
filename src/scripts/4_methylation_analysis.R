#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(R.utils)

where = args[[1]]
nThread = as.numeric(args[[2]])
setwd(where)


# Read all bedgraph files
files = list.files(pattern = "*.cov")
bed_list = lapply(1:length(files), function(x) as.data.frame(fread(files[x], skip = 1, nThread = nThread)))

# Add ID rownames
names(bed_list) = files
for(i in 1:length(bed_list))
{
  rownames(bed_list[[i]]) = paste(bed_list[[i]]$V1, bed_list[[i]]$V2, bed_list[[i]]$V3, sep = "_")
}

# Create a methylation and count matrix across samples
positions = sort(unique(unlist(lapply(bed_list, rownames))))
meth = matrix(nrow = length(positions), ncol = length(files))
coverage = matrix(nrow = length(positions), ncol = length(files))
rownames(meth) = rownames(coverage) = positions
colnames(meth) = colnames(coverage) = files
for(i in 1:length(files))
{
  meth[, i] <- bed_list[[i]][positions,]$V4
  coverage[, i] <- bed_list[[i]][positions,]$V5+bed_list[[i]][positions,]$V6
}

# Export final results
fwrite(data.table(meth, keep.rownames = T), paste(Sys.Date(), "meth.txt", sep = '_'), quote = F,
       row.names = F, col.names = T, sep = '\t', nThread = nThread)

fwrite(data.table(coverage, keep.rownames = T), paste(Sys.Date(), "coverage.txt", sep = '_'), quote = F,
       row.names = F, col.names = T, sep = '\t', nThread = nThread)

setwd("../")
files = list.files(pattern = "*.bedGraph.gz")
bed_list = lapply(1:length(files), function(x) as.data.frame(fread(files[x], skip = 1, nThread = nThread)))




