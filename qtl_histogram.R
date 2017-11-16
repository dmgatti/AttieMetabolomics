################################################################################
# Given the harvested QTL, make a histogram of the top QTL.
# Daniel Gatti
# dan.gatti@jax.org
# July 20, 2017
################################################################################
options(stringsAsFactors = F)
library(tidyverse)

# Arguments:
# input.file: full path to the *.csv qtl summary file.
# output.file: full path to the output figure file as a PNG.
args = commandArgs(trailingOnly = TRUE)
input.file  = args[1]
output.file = args[2]

print(paste("INPUT.FILE=", input.file))
print(paste("OUPUT.FILE=", output.file))

# Load in the data.
data = read.csv(input.file)
data$chr = factor(data$chr, levels = c(1:19, "X"))

# Load in markers.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM.rds")
markers[,2] = factor(markers[,2], levels = c(1:19, "X"))

map = split(markers[,3], markers[,2])
chrlen = sapply(map, length)
chrlen.mb = sapply(map, max)
chrsum = cumsum(chrlen)
chrsum.mb = cumsum(c(0, chrlen.mb[-length(chrlen.mb)]))
names(chrsum.mb) = names(chrlen)
chrmid = c(1, chrsum[-length(chrsum)]) + chrlen / 2
names(chrmid) = names(chrlen)
col = rep(c("grey50", "black"), 10)

# Add genome Mb to data.
data = data.frame(data, gmb = data$pos + chrsum.mb[data$chr])

png(output.file, width = 1200, height = 800, res = 128)

ggplot(data, aes(x = pos)) +
  geom_histogram(binwidth = 2) +
  facet_grid(~chr) + 
  theme(panel.spacing.x = unit(0.1, "lines"))

dev.off()
