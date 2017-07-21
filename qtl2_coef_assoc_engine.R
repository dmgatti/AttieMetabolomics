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
library(qtl2convert)
library(RSQLite)
library(dplyr)

###############################
# Get command line arguments. #
###############################
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
chunk_size   = as.numeric(args[5])
chunk_number = as.numeric(args[6])

data.file  = "/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/qtl2_input/attie_liver_lipids_qtl2_input.Rdata"
qtl.file   = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv"
output.dir = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/"
fig.dir    = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/Liver/lipids_norm_jax/"
chunk_size = 10
chunk_number = 1

################################
# Load in the qtl2 input data. #
################################
load(data.file)

# Load in the QTL summary table.
qtl.summary = read.csv(qtl.file)

stopifnot(c("pheno", "pheno.descr", "genoprobs", "K", "map") %in% ls())

######################
# Set up covariates. #
######################
covar.names = strsplit(pheno.descr$covar_list, ":")
covar.names = unique(unlist(covar.names))
# There may be NA's in the covar.names.
covar.names = covar.names[!is.na(covar.names)]
stopifnot(covar.names %in% colnames(pheno))

f = as.formula(paste("~", paste(covar.names, collapse = "+")))
addcovar = model.matrix(f, data = pheno)[,-1]

#######################################################
# Split out phenotypes and convert to numeric matrix. #
#######################################################

covar = pheno[,!pheno.descr$is_pheno]
pheno = as.matrix(pheno[,pheno.descr$is_pheno])

#########################################
# Calculate the phenotype range to run. #
#########################################

qtl.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
if(qtl.rng[length(qtl.rng)] > nrow(qtl.summary)) {
  qtl.rng = qtl.rng[1]:nrow(qtl.summary)
}

###########################################################################
# Calculate the founder allele effects and association mapping and create #
# plots for each.                                                         #
###########################################################################

# For each analyte, select the peak with the highest LOD, plot the coefs
# on that chromosome and perform association mapping in a +/- 2 MB window
# around the peak.
for(i in qtl.rng) {

  # Get the QTL info for this peak.
  qtl.info = qtl.summary[i,]

  coef = scan1_coef(genoprobs = genoprobs, pheno = pheno[,qtl.info$analyte],
           kinship = K, addcovar = addcovar)

} # for(i)


