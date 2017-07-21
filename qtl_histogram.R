################################################################################
# Make a histogram of the top QTL.
# Daniel Gatti
# dan.gatti@jax.org
# July 20, 2017
################################################################################
options(stringsAsFactors = F)
library(tidyverse)

# Arguments:
# input.file: full path to the *.rds qtl summary file.
# output.file: full path to the output figure file as a PNG.
# thr: LOD threshold to use when selecting QTL peaks.
args = commandArgs(trailingOnly = TRUE)
input.file  = args[1]
output.file = args[2]
thr = as.numeric(args[3])

print(paste("INPUT.FILE=", input.file))
print(paste("OUPUT.FILE=", output.file))
print(paste("THR=", thr))

# Load in the data.
data = readRDS(input.file)

load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[rownames(data), 1:4]
markers[,2] = factor(markers[,2], levels = c(1:19, "X"))

map = split(markers[,3], markers[,2])
chrlen = sapply(map, length)
chrlen = chrlen[order(as.numeric(names(chrlen)))]
chrsum = cumsum(chrlen)
chrmid = c(1, chrsum[-length(chrsum)]) + chrlen / 2
names(chrmid) = names(chrlen)
col = rep(c("grey50", "black"), 10)

png(output.file, width = 1000, height = 800, res = 128)

rs = rowSums(data > thr)
ylim = range(rs)

print(paste("YLIM = ", ylim))

rs.df = data.frame(pos = 1:length(rs), rs = rs)
rs.df = split(rs.df, markers[,2])

par(plt = c(0.08, 0.99, 0.12, 0.95))
plot(-1, -1, type = "l", las = 1, col = col, xlim = c(0, nrow(markers)), 
     xaxt = "n", xlab = "", xaxs = "i", ylim = c(0, ylim[2]),
     ylab = "Number of Analytes")
for(i in 1:length(rs.df)) {
  lines(rs.df[[i]][,1], rs.df[[i]][,2], col = col[i])
}
mtext(side = 1, line = 1, at = chrmid, text = names(chrmid), cex = 1.5)

dev.off()


