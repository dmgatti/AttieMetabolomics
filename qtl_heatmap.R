0###############################################################################
# Given a large matrix with LOD scores and marker map, make a QTL heatmap, 
# clustering the LOD curves, but not the markers.
# Daniel Gatti
# dan.gatti@jax.org
# July 19, 2017
################################################################################
options(stringsAsFactors = F)
library(RColorBrewer)

# Arguments:
# input.file: The *.rds file containing the LOD curves.
# output.file: full path to the figure file to output.
# lod.thr: LOD threshold above which LODs will be truncated to prevent large
#          peaks from dominating the coloring of the heatmap.
args = commandArgs(trailingOnly = TRUE)
input.file  = args[1]
output.file = args[2]
lod.thr = as.numeric(args[3])

print(paste("INPUT.FILE =", input.file))
print(paste("OUTPUT.FILE =", output.file))
print(paste("LOD.THR =", lod.thr))

# Load in the data.
lod = readRDS(input.file)
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[rownames(lod),1:4]

# Get Chr lengths and midpoints.
map = split(markers[,3], markers[,2])
chrlen = sapply(map, length)
chrlen = chrlen[order(as.numeric(names(chrlen)))]
chrsum = cumsum(chrlen)
chrmid = c(1, chrsum[-length(chrsum)]) + chrlen / 2
names(chrmid) = names(chrlen)


# Correlate and cluster the LOD curves.
lod.cor = cor(lod)
cl = hclust(as.dist(1.0 - lod.cor), method = "average")
lod = lod[,cl$order]

# Threshold the maximum LOD to lod.thr so that large QTL don't dominate the heatmap.
lod[lod > lod.thr] = lod.thr
#lod = lod / matrix(colSums(lod), nrow(lod), ncol(lod), byrow = T)

breaks = 0:lod.thr
col = brewer.pal(length(breaks) - 1, "YlOrRd")

png(output.file, width = 6000, height = 4000, res = 128)
par(plt = c(0.1, 0.93, 0.08, 0.95))
image(1:nrow(lod), 1:ncol(lod), lod, breaks = breaks, col = col, axes = F,
      ann = F)
# NOTE: we print every other analyte to make the plot legible.
at = seq(1, ncol(lod), 2)
mtext(side = 2, line = 0.1, at = at, text = colnames(lod)[at], cex = 0.7, las = 1)
mtext(side = 4, line = 0.1, at = at, text = colnames(lod)[at], cex = 0.7, las = 1)
mtext(side = 1, line = 1.5, at = chrmid, text = names(chrmid), cex = 2)
abline(v = chrsum, col = "grey30")
box()
dev.off()
