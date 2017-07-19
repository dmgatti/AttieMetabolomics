################################################################################
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
args = commandArgs(trailingOnly = TRUE)
input.file  = args[1]
output.file = args[2]

print(paste("INPUT.FILE =", input.file))
print(paste("OUTPUT.FILE =", output.file))

# Load in the data.
lod = readRDS(input.file)
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[rownames(lod),1:4]

# Get Chr lengths and midpoints.
map = split(markers[,3], markers[,2])
chrlen = sapply(map, length)


# Correlate and cluster the LOD curves.
lod.cor = cor(lod)
cl = hclust(as.dist(1.0 - lod.cor), method = "average")
lod = lod[,cl$order]





