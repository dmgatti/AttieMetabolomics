################################################################################
# Harvest the QTL for a given data set (i.e. metabolites or lipids).
# Daniel Gatti
# dan.gatti@jas.org
# July 17, 2017
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

input.dir   = args[1]
output.prefix = args[2]
fig.dir     = args[3]

print(paste("INPUT DIR =", input.dir))
print(paste("OUTPUT PREFIX =", output.prefix))
print(paste("FIGURE DIR =", fig.dir))

# Get the QTL files, which may contain more than one result.
qtl.files = dir(path = input.dir, pattern = "_QTL.rds$", full.names = T)

# Read in the marker map and subset to include the markers we used.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
qtl = readRDS(qtl.files[1])
GM_snps = GM_snps[rownames(qtl),1:4]
map = map_df_to_list(map = GM_snps, pos_column = "pos")

qtl.summary = NULL
qtl.mat = NULL

pdf(paste0(fig.dir, "all_QTL.pdf"), width = 10, height = 8)

  for(i in 1:length(qtl.files)) {

    print(paste(i, "of", length(qtl.files)))

    qtl = readRDS(qtl.files[i])

    for(j in 1:ncol(qtl)) {

      plot_scan1(x = qtl, map = map, lodcolumn = j, main = colnames(qtl)[j])

    } # for(j)

    # Harvest and record the maximum QTL.
    max.qtl    = apply(qtl, 2, max)
    wh.max.qtl = apply(qtl, 2, which.max)

    qtl.summary = rbind(qtl.summary, 
                  data.frame(analyte = names(max.qtl),
                             marker = rownames(qtl)[wh.max.qtl],
                             chr = GM_snps$chr[wh.max.qtl],
                             Mb = GM_snps$pos[wh.max.qtl],
                             lod = max.qtl))

    # Add these QTL to the large QTL matrix.
    if(is.null(qtl.mat)) {
      qtl.mat = qtl
    } else {
      qtl.mat = cbind(qtl.mat, qtl)
    } # else
 
  } # for(i)

dev.off()

write.csv(qtl.summary, paste0(output.prefix, "_qtl_summary.csv"), row.names = F,
          quote = F)
saveRDS(qtl.mat, file = paste0(output.prefix, "_all_qtl.rds"))

