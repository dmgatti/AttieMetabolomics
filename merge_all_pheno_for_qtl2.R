################################################################################
# Gather all of the Attie phenotype data and the intersection of the genotype
# data and create one file with everything in it.
# We will use JAX normalized data throughout.
# Use the new data structures from Matt.
# 
# The Rdata file will contain:
# pheno.liver_metabolites: numeric matrix of normalized liver metabolites.
# pheno.liver_lipids: numeric matrix of normalized liver metabolites.
# pheno.cecum_lipids: numeric matrix of normalized liver metabolites.
# pheno.plasma_lipids: numeric matrix of normalized liver metabolites.
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
library(tidyverse)

h5filename = "/hpcdata/gac/derived/CGD_DO_Genoprobs/GigaMUGA_hap_probs_v4.h5"
pheno.dir  = "/hpcdata/gac/derived/Attie_DO_Metabolomics/data/"
pheno.dict.dir = "/hpcdata/gac/raw/Attie_DO_Metabolomics/formatted_data/phenotype_dictionary/"
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

# Phenotypes.
pheno.liver_metab   = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
pheno.liver_lipid   = readRDS(paste0(pheno.dir, "attie_liver_lipids_normalized.rds"))
pheno.cecum_lipid   = readRDS(paste0(pheno.dir, "attie_cecum_lipids_normalized.rds"))
pheno.plasma_lipid  = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))
pheno.islet_protein = readRDS(paste0(pheno.dir, "attie_islet_proteins_normalized.rds"))

rownames(pheno.liver_metab)   = pheno.liver_metab$Mouse.ID
rownames(pheno.liver_lipid)   = pheno.liver_lipid$Mouse.ID
rownames(pheno.cecum_lipid)   = pheno.cecum_lipid$Mouse.ID
rownames(pheno.plasma_lipid)  = pheno.plasma_lipid$Mouse.ID
rownames(pheno.islet_protein) = pheno.islet_protein$Mouse.ID

# Phenotype dictionaries.
annot.liver_metab = read.delim(paste0(pheno.dict.dir, "liver_metabolites_pheno_dict.txt"))
annot.liver_lipid = read.delim(paste0(pheno.dict.dir, "liver_lipids_pheno_dict.txt"))
annot.cecum_lipid = read.delim(paste0(pheno.dict.dir, "cecum_lipids_pheno_dict.txt"))
annot.plasma_lipid = read.delim(paste0(pheno.dict.dir, "plasma_lipids_pheno_dict.txt"))
annot.islet_protein = read.delim(paste0(pheno.dict.dir, "islet_proteins_pheno_dict.txt"))

# Set the colnames of pheno to be the same as the R_name column in the dictionary.
annot.liver_metab$R_name   = colnames(pheno.liver_metab)
annot.liver_lipid$R_name   = colnames(pheno.liver_lipid)
annot.cecum_lipid$R_name   = colnames(pheno.cecum_lipid)
annot.plasma_lipid$R_name  = colnames(pheno.plasma_lipid)
annot.islet_protein$R_name = colnames(pheno.islet_protein)

# Get the union of samples in the phenotypes.
samples = pheno.liver_metab$Mouse.ID
samples = union(samples, pheno.liver_lipid$Mouse.ID)
samples = union(samples, pheno.cecum_lipid$Mouse.ID)
samples = union(samples, pheno.plasma_lipid$Mouse.ID)
samples = union(samples, pheno.islet_protein$Mouse.ID)
samples = intersect(samples, rownames(probs))
samples = sort(samples)

# 457 samples in union.
length(samples)

# Subset all data to contain the same samples.
probs = probs[samples,,]
pheno.liver_metab = pheno.liver_metab[rownames(pheno.liver_metab) %in% samples,]
pheno.liver_lipid = pheno.liver_lipid[rownames(pheno.liver_lipid) %in% samples,]
pheno.cecum_lipid = pheno.cecum_lipid[rownames(pheno.cecum_lipid) %in% samples,]
pheno.plasma_lipid = pheno.plasma_lipid[rownames(pheno.plasma_lipid) %in% samples,]
pheno.islet_protein = pheno.islet_protein[rownames(pheno.islet_protein) %in% samples,]

stopifnot(rownames(pheno.liver_metab) %in% rownames(probs))
stopifnot(rownames(pheno.liver_lipid) %in% rownames(probs))
stopifnot(rownames(pheno.cecum_lipid) %in% rownames(probs))
stopifnot(rownames(pheno.plasma_lipid) %in% rownames(probs))
stopifnot(rownames(pheno.islet_protein) %in% rownames(probs))

# Make covariates (sex, gen & batch).
covar.liver_metab = model.matrix(~sex + wave + Batch, data = pheno.liver_metab)[,-1]
covar.liver_lipid = model.matrix(~sex + wave + Batch, data = pheno.liver_lipid)[,-1]
covar.cecum_lipid = model.matrix(~sex + wave + Batch, data = pheno.cecum_lipid)[,-1]
covar.plasma_lipid = model.matrix(~sex + wave + Batch, data = pheno.plasma_lipid)[,-1]
covar.islet_protein = model.matrix(~sex + wave + Batch, data = pheno.islet_protein)[,-1]

# Make covar_factors (sex, wave, Batch)
covar_factors.liver_metab = data.frame(column_name = c("sex", "wave", "batch"), 
                                       display_name = c("Sex", "Wave", "batch"))
covar_factors.liver_lipid = data.frame(column_name = c("sex", "wave", "batch"), 
                                       display_name = c("Sex", "Wave", "batch"))
covar_factors.cecum_lipid = data.frame(column_name = c("sex", "wave", "batch"), 
                                       display_name = c("Sex", "Wave", "batch"))
covar_factors.plasma_lipid = data.frame(column_name = c("sex", "wave", "batch"), 
                                        display_name = c("Sex", "Wave", "batch"))
covar_factors.islet_protein = data.frame(column_name = c("sex", "wave", "batch"), 
                                         display_name = c("Sex", "Wave", "batch"))

# Make global sample annotation for each data set.
colnames(pheno.liver_metab)   = sub("Batch", "batch", colnames(pheno.liver_metab))
colnames(pheno.liver_metab)   = sub("Prep.Date", "prep_date", colnames(pheno.liver_metab))
colnames(pheno.liver_metab)   = sub("Que.Number", "que_number", colnames(pheno.liver_metab))
colnames(pheno.liver_lipid)   = sub("Batch", "batch", colnames(pheno.liver_lipid))
colnames(pheno.liver_lipid)   = sub("Sample Prep Number", "sample_prep", colnames(pheno.liver_lipid))
colnames(pheno.cecum_lipid)   = sub("Batch", "batch", colnames(pheno.cecum_lipid))
colnames(pheno.cecum_lipid)   = sub("Date", "date", colnames(pheno.cecum_lipid))
colnames(pheno.cecum_lipid)   = sub("Sample_Prep_Order", "sample_prep", colnames(pheno.cecum_lipid))
colnames(pheno.plasma_lipid)  = sub("Batch", "batch", colnames(pheno.plasma_lipid))
colnames(pheno.plasma_lipid)  = sub("Date", "date", colnames(pheno.plasma_lipid))
colnames(pheno.islet_protein) = sub("Batch", "batch", colnames(pheno.islet_protein))
colnames(pheno.islet_protein) = sub("Sample_Prep_Number", "sample_prep", colnames(pheno.islet_protein))
colnames(pheno.islet_protein) = sub("Date", "date", colnames(pheno.islet_protein))

samples.liver_metab   = pheno.liver_metab[,1:13]
samples.liver_lipid   = pheno.liver_lipid[,1:13]
samples.cecum_lipid   = pheno.cecum_lipid[,1:13]
samples.plasma_lipid  = pheno.plasma_lipid[,1:13]
samples.islet_protein = pheno.islet_protein[,1:13]
rownames(samples.liver_metab)   = samples.liver_metab[,1]
rownames(samples.liver_lipid)   = samples.liver_lipid[,1]
rownames(samples.cecum_lipid)   = samples.cecum_lipid[,1]
rownames(samples.plasma_lipid)  = samples.plasma_lipid[,1]
rownames(samples.islet_protein) = samples.islet_protein[,1]
colnames(samples.liver_metab)[1]   = "id"
colnames(samples.liver_lipid)[1]   = "id"
colnames(samples.cecum_lipid)[1]   = "id"
colnames(samples.plasma_lipid)[1]  = "id"
colnames(samples.islet_protein)[1] = "id"

# Subset the phenotypes to only include numeric data and convert them to matrices.
pheno.liver_metab   = as.matrix(pheno.liver_metab[,-(1:13)])
pheno.liver_lipid   = as.matrix(pheno.liver_lipid[,-(1:13)])
pheno.cecum_lipid   = as.matrix(pheno.cecum_lipid[,-(1:13)])
pheno.plasma_lipid  = as.matrix(pheno.plasma_lipid[,-(1:13)])
pheno.islet_protein = as.matrix(pheno.islet_protein[,-(1:13)])

# RankZ transform each phenotype.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

pheno.liver_metab   = apply(pheno.liver_metab, 2, rankZ)
pheno.liver_lipid   = apply(pheno.liver_lipid, 2, rankZ)
pheno.cecum_lipid   = apply(pheno.cecum_lipid, 2, rankZ)
pheno.plasma_lipid  = apply(pheno.plasma_lipid, 2, rankZ)
pheno.islet_protein = apply(pheno.islet_protein, 2, rankZ)

# Get the markers.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
markers = GM_snps[,1:4]

# Subset markers to match probs.
markers = markers[dimnames(probs)[[3]],]
colnames(markers) = c("marker", "chrom", "pos", "cM")

# Convert probs to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = markers, chr_column = "chrom",
            pos_column = "pos")
map = map_df_to_list(markers, chr_column = "chrom", pos_column = "pos")
rm(probs)

# Create LOCO kinship matrices.
K = calc_kinship(probs = genoprobs, type = "loco", cores = 4)

# Load in the pre-computed, harveted LOD peaks.
lod.liver_metab = read.csv("/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_qtl_summary_thresh_6.csv")
lod.liver_lipid = read.csv("/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv")
lod.cecum_lipid = read.csv("/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/lipids_norm_jax/cecum_lipids_jax_norm_qtl_summary_thresh_6.csv")
lod.plasma_lipid = read.csv("/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm_qtl_summary_thresh_6.csv")
lod.islet_protein = read.csv("/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Islet/proteins/islet_proteins_jax_norm_qtl_summary_thresh_6.csv")

colnames(lod.liver_metab)[1]   = "annot.id"
colnames(lod.liver_metab)[2]   = "marker.id"
colnames(lod.liver_lipid)[1]   = "annot.id"
colnames(lod.liver_lipid)[2]   = "marker.id"
colnames(lod.cecum_lipid)[1]   = "annot.id"
colnames(lod.cecum_lipid)[2]   = "marker.id"
colnames(lod.plasma_lipid)[1]  = "annot.id"
colnames(lod.plasma_lipid)[2]  = "marker.id"
colnames(lod.islet_protein)[1] = "annot.id"
colnames(lod.islet_protein)[2] = "marker.id"

# Change the upper case Batch to batch in the phenotype dictionaries.
annot.liver_metab$use_covar   = sub("Batch", "batch", annot.liver_metab$use_covar)
annot.liver_lipid$use_covar   = sub("Batch", "batch", annot.liver_lipid$use_covar)
annot.cecum_lipid$use_covar   = sub("Batch", "batch", annot.cecum_lipid$use_covar)
annot.plasma_lipid$use_covar  = sub("Batch", "batch", annot.plasma_lipid$use_covar)
annot.islet_protein$use_covar = sub("Batch", "batch", annot.islet_protein$use_covar)

# Final QC checks.
stopifnot(rownames(pheno.liver_metab) == rownames(covar.liver_metab))
stopifnot(rownames(pheno.liver_lipid) == rownames(covar.liver_lipid))
stopifnot(rownames(pheno.cecum_lipid) == rownames(covar.cecum_lipid))
stopifnot(rownames(pheno.plasma_lipid) == rownames(covar.plasma_lipid))
stopifnot(rownames(pheno.islet_protein) == rownames(covar.islet_protein))
stopifnot(rownames(pheno.liver_metab) %in% rownames(genoprobs[[1]]))
for(chr in 1:length(K)) {
  stopifnot(rownames(pheno.liver_metab) %in% rownames(K[[chr]]))
  stopifnot(rownames(map)[[chr]] %in% dimnames(genoprobs[[chr]])[[3]])
  stopifnot(rownames(K[[chr]]) == rownames(genoprobs[[chr]]))
} # for(chr)

dim(pheno.liver_metab)
dim(pheno.liver_lipid)
dim(pheno.cecum_lipid)
dim(pheno.plasma_lipid)
dim(pheno.islet_protein)

# Build data sets.
dataset.liver.metab = list(annots    = annot.liver_metab, 
                           covar     = covar.liver_metab,
                           covar.factors = covar_factors.liver_metab,
                           datatype  = "phenotype", 
                           display.name  = "Liver Metabolites",
                           lod.peaks = lod.liver_metab, 
                           pheno     = pheno.liver_metab,
                           samples   = samples.liver_metab)

dataset.liver.lipid = list(annots    = annot.liver_lipid, 
                           covar     = covar.liver_lipid,
                           covar.factors = covar_factors.liver_lipid,
                           datatype  = "phenotype", 
                           display.name  = "Liver Lipids",
                           lod.peaks = lod.liver_lipid, 
                           pheno     = pheno.liver_lipid,
                           samples   = samples.liver_lipid)

dataset.cecum.lipid = list(annots    = annot.cecum_lipid, 
                           covar     = covar.cecum_lipid,
                           covar.factors = covar_factors.cecum_lipid,
                           datatype  = "phenotype", 
                           display.name  = "Cecum Lipids",
                           lod.peaks = lod.cecum_lipid, 
                           pheno     = pheno.cecum_lipid,
                           samples   = samples.cecum_lipid)

dataset.plasma.lipid = list(annots    = annot.plasma_lipid, 
                            covar     = covar.plasma_lipid,
                            covar.factors = covar_factors.plasma_lipid,
                            datatype  = "phenotype", 
                            display.name  = "Plasma Lipids",
                            lod.peaks = lod.plasma_lipid, 
                            pheno     = pheno.plasma_lipid,
                            samples   = samples.plasma_lipid)
 
dataset.islet.protein = list(annots    = annot.islet_protein, 
                             covar     = covar.islet_protein,
                             covar.factors = covar_factors.islet_protein,
                             datatype  = "phenotype", 
                             display.name  = "Islet Proteins",
                             lod.peaks = lod.islet_protein, 
                             pheno     = pheno.islet_protein,
                             samples   = samples.islet_protein)

genome.build = "GRCm38"


# Save the data.
save(dataset.liver.metab, 
     dataset.liver.lipid,
     dataset.cecum.lipid,
     dataset.plasma.lipid,
     dataset.islet.protein,
     genome.build, 
     genoprobs, 
     K,
     map,
     markers,
     file = paste0(output.dir, "attie_all_mass_spec_qtl2_input_v3.Rdata"))


# Run a test scan.
qtl = scan1(genoprobs = genoprobs, pheno = pheno.islet_protein[,1,drop = FALSE],
            kinship = K, addcovar = covar.islet_protein, cores = 4)

plot(qtl, map)


