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
library(AnnotationHub)

source.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/scripts/"

source(paste0(source.dir, "assoc_mapping.R"))

###############################
# Get command line arguments. #
###############################
# Arguments:
# data.file: full path to the qtl2 input data file. *.Rdata file.
# sum.file:  full path to the QTL summary file.
# qtl.file:  full path to the QTL LOD file.
# output.dir: full path to the output directory for the data files.
# fig.dir: full path to the output directory for the data files.
# chunk.size: number of QTL summary rows to map in this run.
# chunk.number: Index of th current chunk of QTL summary rows to run.
args = commandArgs(trailingOnly = TRUE)
data.file  = args[1]
sum.file   = args[2]
qtl.file   = args[3]
output.dir = args[4]
fig.dir    = args[5]
chunk_size   = as.numeric(args[6])
chunk_number = as.numeric(args[7])

#data.file  = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_liver_lipids_qtl2_input.Rdata"
#qtl.file   = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_all_qtl.rds"
#sum.file   = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv"
#output.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/"
#fig.dir    = "/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_lipids_norm_jax/"
#chunk_size = 10
#chunk_number = 1

################################
# Load in the qtl2 input data. #
################################
load(data.file)

# Load in the QTL LOD file.
qtl = readRDS(qtl.file)

# Load in the QTL summary table.
qtl.summary = read.csv(sum.file)

stopifnot(c("pheno", "pheno.descr", "genoprobs", "K", "map") %in% ls())

#############################
# Load in Gigamuga markers. #
#############################
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
marker.names = unlist(lapply(genoprobs, function(z) { dimnames(z)[[3]] }))
markers = GM_snps[marker.names, 1:4]
map = map_df_to_list(markers,  pos_column = "pos")

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

##########################
# Get the Ensembl genes. #
##########################
hub = AnnotationHub()
hub = AnnotationHub::query(hub, c("gtf", "mus musculus", "ensembl"))
ensembl = hub[["AH51040"]]

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
  analyte = qtl.info$analyte
  chr = qtl.info$chr
  pos = qtl.info$pos

  # Founder allele effects.
  coef = scan1coef(genoprobs = genoprobs[,chr], pheno = pheno[,analyte,drop = F],
                   kinship = K[[chr]], addcovar = addcovar)
  saveRDS(coef, file = paste0(output.dir, analyte, "_coef_chr", chr, ".rds"))

  png(paste0(fig.dir, analyte, "_coef_chr", chr, ".png"), width = 1000, height = 1000, res = 128)
  plot_coefCC(x = coef, map = map, scan1_output = qtl[,analyte, drop = FALSE],
              main = analyte)
  dev.off()

  # Association mapping.
  start = max(1, pos - 2)
  end = pos + 2
  assoc = assoc_map(chr = chr, start = start, end = end, probs = genoprobs[,chr], 
                pheno = pheno[, analyte, drop = FALSE], K = K[[chr]], addcovar = addcovar, 
                intcovar = NULL, map = map,  
                db.file = "/data/gac/resource/CCsnps/ccfoundersnps.sqlite", ncores = 4)
  save(assoc, file = paste0(output.dir, analyte, "_assoc_chr", chr, ".Rdata"))

  png(paste0(fig.dir, analyte, "_assoc_chr", chr, ".png"), width = 3000, height = 2000, 
             res = 128)
  assoc_plot(assoc = assoc[[1]], snpinfo = assoc[[2]], map = map, chr = chr, 
             start = start, end = end)
  dev.off()

} # for(i)



