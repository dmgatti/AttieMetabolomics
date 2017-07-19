################################################################################
# Utility functions for mapping Attie mass spec data.
# Daniel Gatti
# dan.gatti@jax.org
# June 17, 2017
################################################################################
options(stringsAsFactors = F)
library(qtl2)   # Loads qtl2geno, qtl2scan & qtl2plot.
library(qtl2convert)
library(RSQLite)
library(dplyr)

# Command line arguments.
# 1: full path to input file name.
# 2: output dir and file prefix.
# 3: chunk_size: integer that is the chunk size.
# 4: chunk_number: integer that is the chunk number to run.
# 5: max_col: integer that is the maximum number of columns in the phenotype file.
args = commandArgs(trailingOnly = TRUE)

input.file = args[1]
output.prefix = args[2]
chunk_size = as.numeric(args[3])
chunk_number = as.numeric(args[4])
max_col = as.numeric(args[5])

print(paste("INPUT:", input.file))
print(paste("OUTPUT:", output.prefix))
print(paste("CHUNK_SIZE:", chunk_size))
print(paste("CHUNK_NUMBER:", chunk_number))
print(paste("MAX_COL:", max_col))

#####################
# Load in the data. #
#####################
load(input.file)

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

pheno.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
if(pheno.rng[length(pheno.rng)] > max_col) {
  pheno.rng = pheno.rng[1]:max_col
}

#########
# Scan1 #
#########
qtl = scan1(genoprobs = genoprobs, pheno = pheno[,pheno.rng,drop = F], kinship = K, 
            addcovar = addcovar, cores = 1)
saveRDS(qtl, file = paste0(output.prefix, "chunk_", chunk_number, "_QTL.rds"))

