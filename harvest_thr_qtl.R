################################################################################
# Harvest the QTL for a given data set.
# Gather together the chunked Q TL files into one LDO file and harvest the
# maximum LOD for each analyte.
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
# 3: LOD threshold.
args = commandArgs(trailingOnly = TRUE)

input.file    = args[1]
output.prefix = args[2]
thr           = as.numeric(args[3])

print(paste("INPUT FILE =", input.file))
print(paste("OUTPUT PREFIX =", output.prefix))
print(paste("THR =", thr))

# Read in the aggregated QTL file.
qtl = readRDS(input.file)

# Read in the marker map and subset to include the markers we used.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[rownames(qtl),1:4]
markers[,2] = factor(markers[,2], levels = c(1:19, "X"))
map = split(markers, markers$chr)

qtl = data.frame(qtl)
rownames(qtl) = markers$marker

qtl= data.frame(marker = markers$marker, chr = markers$chr,
                pos = markers$pos, qtl)
qtl.summary = qtl %>% gather(analyte, lod, -marker, -chr, -pos) %>%
                group_by(analyte, chr) %>%
                filter(lod > thr) %>%
                top_n(1, lod)
qtl.summary = as.data.frame(qtl.summary[,c(4, 1:3, 5)])

write.csv(qtl.summary, paste0(output.prefix, "_qtl_summary_thresh_", thr,".csv"), 
          row.names = F, quote = F)

