################################################################################
# Gather all of the Attie phenotype data and the intersection of the genotype
# data and create one file with everything in it.
# We will use JAX normalize data throughout.
# 
# The Rdata file will contain:
# pheno.liver.metabolites: numeric matrix of normalized liver metabolites.
# pheno.liver.lipids: numeric matrix of normalized liver metabolites.
# pheno.cecum.lipids: numeric matrix of normalized liver metabolites.
# pheno.plasma.lipids: numeric matrix of normalized liver metabolites.
# pheno.islet_proteins: numeric matrix of normalized liver metabolites.
# genoprobs: list of 3D arrays. qtl2-style genoprobs.
# K: list of numeric matrices. LOCO kinship matrices.
# addcovar: matrix of additive covariates.
# markers: data.frame of markers.
# map: list of numeric vectors. qtl2-style marker map.
#
# Daniel Gatti
# dan.gatti@jax.org
# Sept. 14, 2017
################################################################################
options(stringsAsFactors = F)

library(qtl2)
library(qtl2convert)
library(rhdf5)

h5filename = "/hpcdata/gac/derived/CGD_DO_Genoprobs/GigaMUGA_hap_probs_v4.h5"
pheno.dir  = "/hpcdata/gac/derived/Attie_DO_Metabolomics/data/"
output.dir = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/"

# Get genoprobs.
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Load in the phenotypes.
# Liver metabolites.
# Liver lipids.
# Cecum lipids.
# Plasma lipids.
# Islet proteins.

pheno.liver.metab = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
pheno.liver.lipid = readRDS(paste0(pheno.dir, "attie_liver_lipids_normalized.rds"))
pheno.cecum.lipid = readRDS(paste0(pheno.dir, "attie_cecum_lipids_normalized.rds"))
pheno.plasma.lipid = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))
pheno.islet.protein = readRDS(paste0(pheno.dir, "attie_islet_proteins_normalized.rds"))

rownames(pheno.liver.metab) = pheno.liver.metab$Mouse.ID
rownames(pheno.liver.lipid) = pheno.liver.lipid$Mouse.ID
rownames(pheno.cecum.lipid) = pheno.cecum.lipid$Mouse.ID
rownames(pheno.plasma.lipid) = pheno.plasma.lipid$Mouse.ID
rownames(pheno.islet.protein) = pheno.islet.protein$Mouse.ID

# Get the common subset of samples.
samples = rownames(probs)
samples = intersect(samples, pheno.liver.metab$Mouse.ID)
samples = intersect(samples, pheno.liver.lipid$Mouse.ID)
samples = intersect(samples, pheno.cecum.lipid$Mouse.ID)
samples = intersect(samples, pheno.plasma.lipid$Mouse.ID)
samples = intersect(samples, pheno.islet.protein$Mouse.ID)
samples = sort(samples)

# 366 samples
length(samples)

# Subset all data to contain the same samples.
probs = probs[samples,,]
pheno.liver.metab = pheno.liver.metab[samples,]
pheno.liver.lipid = pheno.liver.lipid[samples,]
pheno.cecum.lipid = pheno.cecum.lipid[samples,]
pheno.plasma.lipid = pheno.plasma.lipid[samples,]
pheno.islet.protein = pheno.islet.protein[samples,]

stopifnot(rownames(probs) == rownames(pheno.liver.metab))
stopifnot(rownames(probs) == rownames(pheno.liver.lipid))
stopifnot(rownames(probs) == rownames(pheno.cecum.lipid))
stopifnot(rownames(probs) == rownames(pheno.plasma.lipid))
stopifnot(rownames(probs) == rownames(pheno.islet.protein))

# Make covariates (sex, gen & batch).
covar.liver.metab = model.matrix(~sex + wave + Batch, data = pheno.liver.metab)[,-1]
covar.liver.lipid = model.matrix(~sex + wave + Batch, data = pheno.liver.lipid)[,-1]
covar.cecum.lipid = model.matrix(~sex + wave + Batch, data = pheno.cecum.lipid)[,-1]
covar.plasma.lipid = model.matrix(~sex + wave + Batch, data = pheno.plasma.lipid)[,-1]
covar.islet.protien = model.matrix(~sex + wave + Batch, data = pheno.islet.protien)[,-1]

# Subset the phenotypes to only include numeric data and convert them to matrices.
pheno.liver.metab = as.matrix(pheno.liver.metab[,-(1:13)])
pheno.liver.lipid = as.matrix(pheno.liver.lipid[,-(1:13)])
pheno.cecum.lipid = as.matrix(pheno.cecum.lipid[,-(1:13)])
pheno.plasma.lipid = as.matrix(pheno.plasma.lipid[,-(1:13)])
pheno.islet.protein = as.matrix(pheno.islet.protein[,-(1:13)])

# RankZ transform each phenotype.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

pheno.liver.metab = apply(pheno.liver.metab, 2, rankZ)
pheno.liver.lipid = apply(pheno.liver.lipid, 2, rankZ)
pheno.cecum.lipid = apply(pheno.cecum.lipid, 2, rankZ)
pheno.plasma.lipid = apply(pheno.plasma.lipid, 2, rankZ)
pheno.islet.protein = apply(pheno.islet.protein, 2, rankZ)

# Get the markers.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[,1:4]

# Subset markers to match probs.
markers = markers[dimnames(probs)[[3]],]

# Convert probs to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = markers, pos_column = "pos")
map = map_df_to_list(markers, pos_column = "pos")
rm(probs)

# Create LOCO kinship matrices.
K = calc_kinship(probs = genoprobs, type = "loco", cores = 4)

# Final QC checks.
stopifnot(rownames(pheno.liver.metab) == rownames(covar.liver.metab))
stopifnot(rownames(pheno.liver.lipid) == rownames(covar.liver.lipid))
stopifnot(rownames(pheno.cecum.lipid) == rownames(covar.cecum.lipid))
stopifnot(rownames(pheno.plasma.lipid) == rownames(covar.plasma.lipid))
stopifnot(rownames(pheno.islet.protein) == rownames(covar.islet.protein))
stopifnot(rownames(pheno.liver.metab) == rownames(genoprobs[[1]]))
for(chr in 1:length(K)) {
  stopifnot(rownames(pheno.liver.metab) == rownames(K[[chr]]))
  stopifnot(rownames(map)[[chr]] == dimnames(genoprobs[[chr]])[[3]])
} # for(chr)

dim(pheno.liver.metab)
dim(pheno.liver.lipid)
dim(pheno.cecum.lipid)
dim(pheno.plasma.lipid)
dim(pheno.islet.protein)

# Save the data.
save(pheno.liver.metab, pheno.liver.lipid, pheno.cecum.lipid, pheno.plasma.lipid,
     covar.liver.metab, covar.liver.lipid, covar.cecum.lipid, covar.plasma.lipid,
     genoprobs, K, markers, map, file = paste0(output.dir, "attie_all_mass_spec_qtl2_input.Rdata"))

