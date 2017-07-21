################################################################################
# Gather together the Attie liver metabolite, lipid, covariate and genoprobs
# data. Create objects taht can be fed into qtl2 for mapping including:
# pheno: data.frame containing covariates and phenotypes.
# pheno.descr: data.frame with pheotype descriptions.
# genoprobs: qtl2 style genoprobs object, one per chromosome.
# K: list containing kinship matrices, one per chromosome.
# map: list of marker positions, one per chromosome.
#
# We will create two Rdata files: one each for metabolites and lipids.
# 
# Daniel Gatti
# dan.gatti@jax.org
# July 17, 2017
################################################################################
options(stringsAsFactors = F)
library(qtl2)
library(qtl2convert)
library(rhdf5)
library(tidyverse)

# Set up input directories.
pheno.dir = "/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/data/"
probs.dir = "/hpcdata/gac/derived/CGD_DO_Genoprobs/"
uwisc.dir = "/hpcdata/gac/raw/Attie_DO_Liver_Metabolomics/formatted_data/"
out.dir   = "/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/qtl2_input/"

#########################################################################
# Liver Metabolites: JAX normalization: sex, generation and batch covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("metabolites", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(probs)))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
probs    = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_metabolites_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#################################################################
# Liver Metabolites: JAX normalization, sex and generation covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("metabolites", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(probs)))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
probs    = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_metabolites_qtl2_input_sex_gen.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)



######################################################################
# Liver Lipids: JAX normalization, sex, generation & batch covariates.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_liver_lipids_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("lipids", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_lipids_zscore_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(probs)))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
probs    = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_lipids_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#######################################################################
# Plasma Lipids: JAX normalization, sex, generation & batch covariates.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("lipids", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_plasma_lipids_zscore_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(probs)))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
probs    = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_plasma_lipids_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


##########################################################################
# Metabolites: U. Wisc. normalization, sex, generation & batch covariates.

# Load in the phenotypes.
pheno = read.delim(paste0(uwisc.dir, "14June2017_DOLiverMetabolites_NORM.txt"))
rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
covar = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))

pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
rownames(pheno) = pheno$Mouse.ID
rm(covar)

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("metabolites", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples = sort(intersect(rownames(pheno), rownames(probs)))
if(length(samples) == 0) stop("NO SAMPLES IN COMMON!!")
pheno = pheno[samples,]
probs = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_metabolites_UWisc_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


##########################################################################
# Metabolites: U. Wisc. normalization, sex & generation covariates.

# Load in the phenotypes.
pheno = read.delim(paste0(uwisc.dir, "14June2017_DOLiverMetabolites_NORM.txt"))
rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
covar = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))

pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
rownames(pheno) = pheno$Mouse.ID
rm(covar)

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("metabolites", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples = sort(intersect(rownames(pheno), rownames(probs)))
if(length(samples) == 0) stop("NO SAMPLES IN COMMON!!")
pheno = pheno[samples,]
probs = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_metabolites_UWisc_qtl2_input_sex_gen.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#######################################
# Liver Lipids: U. Wisc. normalization.

# Load in the phenotypes.
pheno = read.delim(paste0(uwisc.dir, "15June2017_DOLiverLipidomics.txt"))
rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
covar = readRDS(paste0(pheno.dir, "attie_liver_lipids_normalized.rds"))

pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
rownames(pheno) = pheno$Mouse.ID
rm(covar)

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("lipids", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_lipids_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples = sort(intersect(rownames(pheno), rownames(probs)))
if(length(samples) == 0) stop("NO SAMPLES IN COMMON!!")
pheno = pheno[samples,]
probs = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_lipids_qtl2_UWisc_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)



#######################################
# Plasma Lipids: U. Wisc. normalization.

# Load in the phenotypes.
pheno = read.delim(paste0(uwisc.dir, "15June2017_DOPlasmaLipidomics.txt"))
rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
covar = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))

pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
rownames(pheno) = pheno$Mouse.ID
rm(covar)

# Create data dictionary.
num.covar.columns = 13
num.pheno.columns = ncol(pheno) - 13
pheno.descr = data.frame(
              name = colnames(pheno),
              short_name = colnames(pheno),
              description = colnames(pheno),
              category = c(rep(NA, num.covar.columns), rep("lipids", num.pheno.columns)),
              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
              )

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_plasma_lipids_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
h5filename = paste0(probs.dir, "GigaMUGA_hap_probs_v2.h5")
h5ls(h5filename)
grp = h5read(file = h5filename, name = "/ADA")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
rownames(probs) = sub("^AA\\.", "", rownames(probs))
rownames(probs) = sub("(f|m)$", "", rownames(probs))
rownames(probs) = sub("^DO", "DO-", rownames(probs))

# Subset by samples.
samples = sort(intersect(rownames(pheno), rownames(probs)))
if(length(samples) == 0) stop("NO SAMPLES IN COMMON!!")
pheno = pheno[samples,]
probs = probs[samples,,]

print(paste(nrow(pheno), "Samples"))

# Read in the marker data.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))

# Make marker map.
map = map_df_to_list(map = GM_snps[,1:4], pos_column="pos")

# Convert to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = probs, map = GM_snps[,1:4], pos_column="pos")
rm(probs)

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_plasma_lipids_qtl2_UWisc_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)

