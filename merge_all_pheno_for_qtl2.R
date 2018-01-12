################################################################################
# Gather all of the Attie phenotype data and the intersection of the genotype
# data and create one file with everything in it.
# We will use JAX normalized data throughout.
#
# Convert all sample IDs to DOnn format.
#
# Use the new data structures from Matt at:
#  https://github.com/churchill-lab/qtl-viewer/blob/master/docs/QTLViewerDataStructures.md
# 
# The Rdata file will contain:
# 
# genome.build - string specifying the genome build
# genoprobs - the genotype probabilities
# K - the kinship matrix
# map - list of one element per chromosome, with the genomic position of each marker
# markers - marker names and positions

# dataset.liver_metabolites
# dataset.liver_lipids
# dataset.cecum_metabolites
# dataset.cecum_lipids
# dataset.plasma_lipids
# dataset.plasma_lipids
# dataset.islet_proteins
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

genoprobs.file = "/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds"
pheno.dir  = "/hpcdata/gac/derived/Attie_DO_Metabolomics/data/"
pheno.dict.dir = "/hpcdata/gac/raw/Attie_DO_Metabolomics/formatted_data/phenotype_dictionary/"
output.dir = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/"

# Get genoprobs.
genoprobs = readRDS(genoprobs.file)

# Load in the phenotypes.
# Liver metabolites.
# Liver lipids.
# Cecum metabolites.
# Cecum lipids.
# Plasma metabolites.
# Plasma lipids.
# Islet proteins.

# Phenotypes.
pheno.liver_metab   = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
pheno.liver_lipid   = readRDS(paste0(pheno.dir, "attie_liver_lipids_normalized.rds"))
pheno.cecum_metab   = readRDS(paste0(pheno.dir, "attie_cecum_metabolites_normalized.rds"))
pheno.cecum_lipid   = readRDS(paste0(pheno.dir, "attie_cecum_lipids_normalized.rds"))
pheno.plasma_metab  = readRDS(paste0(pheno.dir, "attie_plasma_metabolites_normalized.rds"))
pheno.plasma_lipid  = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))
pheno.islet_protein = readRDS(paste0(pheno.dir, "attie_islet_proteins_normalized.rds"))

colnames(pheno.liver_metab)[1]   = "id"
colnames(pheno.liver_lipid)[1]   = "id"
colnames(pheno.cecum_metab)[1]   = "id"
colnames(pheno.cecum_lipid)[1]   = "id"
colnames(pheno.plasma_metab)[1]  = "id"
colnames(pheno.plasma_lipid)[1]  = "id"
colnames(pheno.islet_protein)[1] = "id"

# Convet sample IDs to format without dashes or periods.
pheno.liver_metab$id   = sub("\\-|\\.", "", pheno.liver_metab$id)
pheno.liver_lipid$id   = sub("\\-|\\.", "", pheno.liver_lipid$id)
pheno.cecum_metab$id   = sub("\\-|\\.", "", pheno.cecum_metab$id)
pheno.cecum_lipid$id   = sub("\\-|\\.", "", pheno.cecum_lipid$id)
pheno.plasma_metab$id  = sub("\\-|\\.", "", pheno.plasma_metab$id)
pheno.plasma_lipid$id  = sub("\\-|\\.", "", pheno.plasma_lipid$id)
pheno.islet_protein$id = sub("\\-|\\.", "", pheno.islet_protein$id)

rownames(pheno.liver_metab)   = pheno.liver_metab$id
rownames(pheno.liver_lipid)   = pheno.liver_lipid$id
rownames(pheno.cecum_metab)   = pheno.cecum_metab$id
rownames(pheno.cecum_lipid)   = pheno.cecum_lipid$id
rownames(pheno.plasma_metab)  = pheno.plasma_metab$id
rownames(pheno.plasma_lipid)  = pheno.plasma_lipid$id
rownames(pheno.islet_protein) = pheno.islet_protein$id

# Phenotype dictionaries.
annot.liver_metab = read.delim(paste0(pheno.dict.dir, "liver_metabolites_pheno_dict.txt"))
annot.liver_lipid = read.delim(paste0(pheno.dict.dir, "liver_lipids_pheno_dict.txt"))
annot.cecum_metab = read.delim(paste0(pheno.dict.dir, "cecum_metabolites_pheno_dict.txt"))
annot.cecum_lipid = read.delim(paste0(pheno.dict.dir, "cecum_lipids_pheno_dict.txt"))
annot.plasma_metab = read.delim(paste0(pheno.dict.dir, "plasma_metabolites_pheno_dict.txt"))
annot.plasma_lipid = read.delim(paste0(pheno.dict.dir, "plasma_lipids_pheno_dict.txt"))
annot.islet_protein = read.delim(paste0(pheno.dict.dir, "islet_proteins_pheno_dict.txt"))

annot.liver_metab[1,1:3]   = "id"
annot.liver_lipid[1,1:3]   = "id"
annot.cecum_metab[1,1:3]   = "id"
annot.cecum_lipid[1,1:3]   = "id"
annot.plasma_metab[1,1:3]  = "id"
annot.plasma_lipid[1,1:3]  = "id"
annot.islet_protein[1,1:3] = "id"

# Change sex, batch and wave column names to be consistent.
colnames(pheno.liver_metab)[grep("sex", colnames(pheno.liver_metab), ignore.case = TRUE)] = "sex"
colnames(pheno.liver_lipid)[grep("sex", colnames(pheno.liver_lipid), ignore.case = TRUE)] = "sex"
colnames(pheno.cecum_metab)[grep("sex", colnames(pheno.cecum_metab), ignore.case = TRUE)] = "sex"
colnames(pheno.cecum_lipid)[grep("sex", colnames(pheno.cecum_lipid), ignore.case = TRUE)] = "sex"
colnames(pheno.plasma_metab)[grep("sex", colnames(pheno.plasma_metab), ignore.case = TRUE)] = "sex"
colnames(pheno.plasma_lipid)[grep("sex", colnames(pheno.plasma_lipid), ignore.case = TRUE)] = "sex"
colnames(pheno.islet_protein)[grep("sex", colnames(pheno.islet_protein), ignore.case = TRUE)] = "sex"

colnames(pheno.liver_metab)[grep("wave", colnames(pheno.liver_metab), ignore.case = TRUE)] = "DOwave"
colnames(pheno.liver_lipid)[grep("wave", colnames(pheno.liver_lipid), ignore.case = TRUE)] = "DOwave"
colnames(pheno.cecum_metab)[grep("wave", colnames(pheno.cecum_metab), ignore.case = TRUE)] = "DOwave"
colnames(pheno.cecum_lipid)[grep("wave", colnames(pheno.cecum_lipid), ignore.case = TRUE)] = "DOwave"
colnames(pheno.plasma_metab)[grep("wave", colnames(pheno.plasma_metab), ignore.case = TRUE)] = "DOwave"
colnames(pheno.plasma_lipid)[grep("wave", colnames(pheno.plasma_lipid), ignore.case = TRUE)] = "DOwave"
colnames(pheno.islet_protein)[grep("wave", colnames(pheno.islet_protein), ignore.case = TRUE)] = "DOwave"

colnames(pheno.liver_metab)[grep("batch", colnames(pheno.liver_metab), ignore.case = TRUE)] = "batch"
colnames(pheno.liver_lipid)[grep("batch", colnames(pheno.liver_lipid), ignore.case = TRUE)] = "batch"
colnames(pheno.cecum_metab)[grep("batch", colnames(pheno.cecum_metab), ignore.case = TRUE)] = "batch"
colnames(pheno.cecum_lipid)[grep("batch", colnames(pheno.cecum_lipid), ignore.case = TRUE)] = "batch"
colnames(pheno.plasma_metab)[grep("batch", colnames(pheno.plasma_metab), ignore.case = TRUE)] = "batch"
colnames(pheno.plasma_lipid)[grep("batch", colnames(pheno.plasma_lipid), ignore.case = TRUE)] = "batch"
colnames(pheno.islet_protein)[grep("batch", colnames(pheno.islet_protein), ignore.case = TRUE)] = "batch"

# Make global sample annotation for each data set.
colnames(pheno.liver_lipid)   = sub("Batch", "batch", colnames(pheno.liver_lipid))
colnames(pheno.liver_lipid)   = sub("Sample Prep Number", "Sample_Prep_Order", colnames(pheno.liver_lipid))
colnames(pheno.liver_lipid)   = sub("Date", "date", colnames(pheno.liver_lipid))
colnames(pheno.cecum_lipid)   = sub("Batch", "batch", colnames(pheno.cecum_lipid))
colnames(pheno.cecum_lipid)   = sub("Date", "date", colnames(pheno.cecum_lipid))
colnames(pheno.plasma_lipid)  = sub("Batch", "batch", colnames(pheno.plasma_lipid))
colnames(pheno.plasma_lipid)  = sub("Date", "date", colnames(pheno.plasma_lipid))
colnames(pheno.plasma_lipid)  = sub("Sample Prep Number", "Sample_Prep_Order", colnames(pheno.plasma_lipid))

# Set the colnames of pheno to be the same as the R_name column in the dictionary.
annot.liver_metab$data_name   = colnames(pheno.liver_metab)
annot.liver_lipid$data_name   = colnames(pheno.liver_lipid)
annot.cecum_metab$data_name   = colnames(pheno.cecum_metab)
annot.cecum_lipid$data_name   = colnames(pheno.cecum_lipid)
annot.plasma_metab$data_name  = colnames(pheno.plasma_metab)
annot.plasma_lipid$data_name  = colnames(pheno.plasma_lipid)
annot.islet_protein$data_name = colnames(pheno.islet_protein)

annot.liver_metab$R_name   = colnames(pheno.liver_metab)
annot.liver_lipid$R_name   = colnames(pheno.liver_lipid)
annot.cecum_metab$R_name   = colnames(pheno.cecum_metab)
annot.cecum_lipid$R_name   = colnames(pheno.cecum_lipid)
annot.plasma_metab$R_name  = colnames(pheno.plasma_metab)
annot.plasma_lipid$R_name  = colnames(pheno.plasma_lipid)
annot.islet_protein$R_name = colnames(pheno.islet_protein)

# Make sure that ll of the dictionaries use the dame covariates.
annot.liver_metab$use_covar   = "sex:DOwave:batch"
annot.liver_lipid$use_covar   = "sex:DOwave:batch"
annot.cecum_metab$use_covar   = "sex:DOwave:batch"
annot.cecum_lipid$use_covar   = "sex:DOwave:batch"
annot.plasma_metab$use_covar  = "sex:DOwave:batch"
annot.plasma_lipid$use_covar  = "sex:DOwave:batch"
annot.islet_protein$use_covar = "sex:DOwave:batch"

# Get the union of samples in the phenotypes.
samples = pheno.liver_metab$id
samples = union(samples, pheno.liver_lipid$id)
samples = union(samples, pheno.cecum_metab$id)
samples = union(samples, pheno.cecum_lipid$id)
samples = union(samples, pheno.plasma_metab$id)
samples = union(samples, pheno.plasma_lipid$id)
samples = union(samples, pheno.islet_protein$id)
samples = intersect(samples, rownames(genoprobs[[1]]))
samples = sort(samples)

# 472 samples in union.
length(samples)

stopifnot(rownames(pheno.liver_metab)   %in% rownames(genoprobs[[1]]))
stopifnot(rownames(pheno.liver_lipid)   %in% rownames(genoprobs[[1]]))
stopifnot(rownames(pheno.cecum_metab)   %in% rownames(genoprobs[[1]]))
stopifnot(rownames(pheno.cecum_lipid)   %in% rownames(genoprobs[[1]]))
stopifnot(rownames(pheno.plasma_metab)  %in% rownames(genoprobs[[1]]))
stopifnot(rownames(pheno.plasma_lipid)  %in% rownames(genoprobs[[1]]))
stopifnot(rownames(pheno.islet_protein) %in% rownames(genoprobs[[1]]))

# Make covariates (sex, wave & batch).
covar.liver_metab   = model.matrix(~sex + DOwave + batch, data = pheno.liver_metab)[,-1]
covar.liver_lipid   = model.matrix(~sex + DOwave + batch, data = pheno.liver_lipid)[,-1]
covar.cecum_metab   = model.matrix(~sex + DOwave + batch, data = pheno.cecum_metab)[,-1]
covar.cecum_lipid   = model.matrix(~sex + DOwave + batch, data = pheno.cecum_lipid)[,-1]
covar.plasma_metab  = model.matrix(~sex + DOwave + batch, data = pheno.plasma_metab)[,-1]
covar.plasma_lipid  = model.matrix(~sex + DOwave + batch, data = pheno.plasma_lipid)[,-1]
covar.islet_protein = model.matrix(~sex + DOwave + batch, data = pheno.islet_protein)[,-1]

# Make covar_factors (sex, wave, batch)
covar_factors.liver_metab = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                       display_name = c("Sex", "Wave", "Batch"))
covar_factors.liver_lipid = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                       display_name = c("Sex", "Wave", "Batch"))
covar_factors.cecum_metab = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                       display_name = c("Sex", "Wave", "Batch"))
covar_factors.cecum_lipid = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                       display_name = c("Sex", "Wave", "Batch"))
covar_factors.plasma_metab = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                        display_name = c("Sex", "Wave", "Batch"))
covar_factors.plasma_lipid = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                        display_name = c("Sex", "Wave", "Batch"))
covar_factors.islet_protein = data.frame(column_name = c("sex", "DOwave", "batch"), 
                                         display_name = c("Sex", "Wave", "Batch"))

samples.liver_metab   = data.frame(pheno.liver_metab[,1:11])
samples.liver_lipid   = data.frame(pheno.liver_lipid[,1:13])
samples.cecum_metab   = data.frame(pheno.cecum_metab[,1:11])
samples.cecum_lipid   = data.frame(pheno.cecum_lipid[,1:13])
samples.plasma_metab  = data.frame(pheno.plasma_metab[,1:11])
samples.plasma_lipid  = data.frame(pheno.plasma_lipid[,1:13])
samples.islet_protein = data.frame(pheno.islet_protein[,1:13])

rownames(samples.liver_metab)   = samples.liver_metab[,1]
rownames(samples.liver_lipid)   = samples.liver_lipid[,1]
rownames(samples.cecum_metab)   = samples.cecum_metab[,1]
rownames(samples.cecum_lipid)   = samples.cecum_lipid[,1]
rownames(samples.plasma_metab)  = samples.plasma_metab[,1]
rownames(samples.plasma_lipid)  = samples.plasma_lipid[,1]
rownames(samples.islet_protein) = samples.islet_protein[,1]

colnames(samples.liver_metab)[1]   = "id"
colnames(samples.liver_lipid)[1]   = "id"
colnames(samples.cecum_metab)[1]   = "id"
colnames(samples.cecum_lipid)[1]   = "id"
colnames(samples.plasma_metab)[1]  = "id"
colnames(samples.plasma_lipid)[1]  = "id"
colnames(samples.islet_protein)[1] = "id"

# RankZ transform each phenotype.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

pheno.liver_metab[,-(1:11)]   = apply(pheno.liver_metab[,-(1:11)],   2, rankZ)
pheno.liver_lipid[,-(1:13)]   = apply(pheno.liver_lipid[,-(1:13)],   2, rankZ)
pheno.cecum_metab[,-(1:11)]   = apply(pheno.cecum_metab[,-(1:11)],   2, rankZ)
pheno.cecum_lipid[,-(1:13)]   = apply(pheno.cecum_lipid[,-(1:13)],   2, rankZ)
pheno.plasma_metab[,-(1:11)]  = apply(pheno.plasma_metab[,-(1:11)],  2, rankZ)
pheno.plasma_lipid[,-(1:13)]  = apply(pheno.plasma_lipid[,-(1:13)],  2, rankZ)
pheno.islet_protein[,-(1:13)] = apply(pheno.islet_protein[,-(1:13)], 2, rankZ)

# Substitute the ENSMUST IDs given by U. Wisc. for the UNIPROT IDs.
u2e = read.delim("/hpcdata/gac/raw/Attie_DO_Metabolomics/raw_data/Uniprot2Ensembl.tab")
spl = strsplit(colnames(pheno.islet_protein), "_")
m = lapply(spl, match, u2e$From)
m = lapply(m, function(z) { z[!is.na(z)] })

mapping = colnames(pheno.islet_protein)

wh = which(sapply(m, length) > 0)

for(i in wh) {

  print(i)
  mapping[i] = paste(unique(u2e$To[m[[i]]]), collapse = "_")

} # for(i)

colnames(pheno.islet_protein) = mapping
annot.islet_protein[,1] = mapping
annot.islet_protein[,2] = mapping
annot.islet_protein[,3] = mapping

# Get the markers.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make sure that all of the markers match between genoprobs and markers.
stopifnot(rownames(markers) == unlist(sapply(genoprobs, function(z) { dimnames(z)[[3]] })))

# Convert markers to qtl2 format.
map = map_df_to_list(markers, chr_column = "chr", pos_column = "pos")

# Create LOCO kinship matrices.
K = calc_kinship(probs = genoprobs, type = "loco", cores = 4)

# Load in the pre-computed, harveted LOD peaks.
base.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/"
lod.liver_metab = read.csv(paste0(base.dir, "Liver/metabolites_norm_jax/liver_metabolites_jax_norm_qtl_summary_thresh_6.csv"))
lod.liver_lipid = read.csv(paste0(base.dir, "Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv"))
lod.cecum_metab = read.csv(paste0(base.dir, "Cecum/metabolites/cecum_metabolites_qtl_summary_thresh_6.csv"))
lod.cecum_lipid = read.csv(paste0(base.dir,  "Cecum/lipids_norm_jax/cecum_lipids_jax_norm_qtl_summary_thresh_6.csv"))
lod.plasma_metab = read.csv(paste0(base.dir, "Plasma/metabolites/plasma_metabolites_qtl_summary_thresh_6.csv"))
lod.plasma_lipid = read.csv(paste0(base.dir, "Plasma/lipids_norm_jax/plasma_lipids_jax_norm_qtl_summary_thresh_6.csv"))
lod.islet_protein = read.csv(paste0(base.dir, "Islet/proteins/islet_proteins_jax_norm_qtl_summary_thresh_6.csv"))

colnames(lod.liver_metab)[1:2]   = c("annot.id", "marker.id")
colnames(lod.liver_lipid)[1:2]   = c("annot.id", "marker.id")
colnames(lod.cecum_metab)[1:2]   = c("annot.id", "marker.id")
colnames(lod.cecum_lipid)[1:2]   = c("annot.id", "marker.id")
colnames(lod.plasma_metab)[1:2]  = c("annot.id", "marker.id")
colnames(lod.plasma_lipid)[1:2]  = c("annot.id", "marker.id")
colnames(lod.islet_protein)[1:2] = c("annot.id", "marker.id")

# Change the upper case Batch to batch in the phenotype dictionaries.
annot.liver_metab$use_covar   = sub("Batch", "batch", annot.liver_metab$use_covar)
annot.liver_lipid$use_covar   = sub("Batch", "batch", annot.liver_lipid$use_covar)
annot.cecum_metab$use_covar   = sub("Batch", "batch", annot.cecum_metab$use_covar)
annot.cecum_lipid$use_covar   = sub("Batch", "batch", annot.cecum_lipid$use_covar)
annot.plasma_metab$use_covar  = sub("Batch", "batch", annot.plasma_metab$use_covar)
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
dim(pheno.cecum_metab)
dim(pheno.cecum_lipid)
dim(pheno.plasma_metab)
dim(pheno.plasma_lipid)
dim(pheno.islet_protein)

source("/hpcdata/gac/projects/Attie_DO_Metabolomics/scripts/qtlDataCheck.R")

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

CheckVariables()
CheckDatasets()


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
qtl = scan1(genoprobs = genoprobs, pheno = pheno.islet_protein[,20,drop = FALSE],
            kinship = K, addcovar = covar.islet_protein, cores = 4)

plot(qtl, map)


