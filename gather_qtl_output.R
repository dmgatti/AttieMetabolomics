################################################################################
# Gather the chunked QTL files for a given data set (i.e. metabolites or lipids).
# Also make a large PDF containing all QTL plots.
# Daniel Gatti
# dan.gatti@jas.org
# July 21, 2017
################################################################################
options(stringsAsFactors = F)
library(qtl2)
library(qtl2convert)
library(tidyverse)

# Command line arguments.
# 1: full path to QTL directory.
# 2: full path to output file prefix. We append .csv and .rds to it.
# 3: full path to figure directory.
args = commandArgs(trailingOnly = TRUE)

input.dir     = args[1]
output.prefix = args[2]
fig.dir       = args[3]

print(paste("INPUT DIR =", input.dir))
print(paste("OUTPUT PREFIX =", output.prefix))
print(paste("FIGURE DIR =", fig.dir))

# Get the QTL files, which may contain more than one result.
qtl.files = dir(path = input.dir, pattern = "_QTL.rds$", full.names = T)

# Read in the marker map and subset to include the markers we used.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM.rds")
map = map_df_to_list(map = markers, pos_column = "pos")

qtl.mat = NULL

pdf(paste0(fig.dir, "all_QTL.pdf"), width = 10, height = 8)

  for(i in 1:length(qtl.files)) {

    print(paste(i, "of", length(qtl.files)))

    qtl = readRDS(qtl.files[i])

    for(j in 1:ncol(qtl)) {

      plot_scan1(x = qtl, map = map, lodcolumn = j, main = colnames(qtl)[j])

    } # for(j)

    # Add these QTL to the large QTL matrix.
    if(is.null(qtl.mat)) {
      qtl.mat = qtl
    } else {
      qtl.mat = cbind(qtl.mat, qtl)
    } # else
 
  } # for(i)

dev.off()

class(qtl.mat) = c("scan1", "matrix")
saveRDS(qtl.mat, file = paste0(output.prefix, "_all_qtl.rds"))

