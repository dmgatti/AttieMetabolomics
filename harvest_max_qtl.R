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
# 1: full path to aggregated QTL LOD file (as *.rds).
# 2: full path to output file prefix. We append .csv and .rds to it.
args = commandArgs(trailingOnly = TRUE)

input.file   = args[1]
output.prefix = args[2]

print(paste("INPUT DIR =", input.file))
print(paste("OUTPUT PREFIX =", output.prefix))

# Load in the LOD file.
qtl = readRDS(input.file)

# Read in the marker map and subset to include the markers we used.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[rownames(qtl),1:4]
markers$chr = factor(markers$chr, levels = c(1:19, "X"))
map = map_df_to_list(map = markers, pos_column = "pos")

qtl.summary = NULL

qtl = data.frame(marker = markers$marker, chr = markers$chr, pos = markers$pos, qtl)

qtl.summary = qtl %>% gather(analyte, lod, -marker, -chr, -pos) %>%
                group_by(analyte) %>%
                top_n(1, lod) %>%
                arrange(desc(lod))

qtl.summary = qtl.summary[,c(4, 1:3, 5)]

write.csv(qtl.summary, paste0(output.prefix, "_qtl_summary.csv"), row.names = F,
          quote = F)


