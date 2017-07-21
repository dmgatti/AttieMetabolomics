################################################################################
# Given the QTL mapping input data and a table of peaks from one of the QTL
# harvest functions, calculate the founder allele effects and association 
# mapping at each peak.
# 
# Daniel Gatti
# dan.gatti@jax.org
# July 20, 2017
################################################################################
options(stringsAsFactors = F)
library(qtl2)

# Arguments:
# data.file: full path to the qtl2 input data file. *.Rdata file.
# qtl.file:  full path to the QTL summary file.
# output.dir: full path to the output directory for the data files.
# fig.dir: full path to the output directory for the data files.
# chunk.size: number of QTL summary rows to map in this run.
# chunk.number: Index of th current chunk of QTL summary rows to run.

args = commandArgs(trailingOnly = TRUE)
data.file  = args[1]
qtl.file   = args[2]
output.dir = args[3]
fig.dir    = args[4]
chunk.size   = as.numeric(args[5])
chunk.number = as.numeric(args[6])

data.file  = "/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/qtl2_input/attie_liver_lipids_qt> _input.Rdata"
qtl.file   = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv"
output.dir = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/"
fig.dir    = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/Liver/lipids_norm_jax/"
chunk.size = 10
chunk.number = 1


# Load in the qtl2 input data.
load(data.file)

# Load in the QTL summary table.
qtl.summary = read.csv(qtl.file)

# For each analyte, select the peak with the highest LOD, plot the coefs
# on that chromosome and perform association mapping in a +/- 2 MB window
# around the peak.
for(i in 1:ncol(qtl)) {




} # for(i)


