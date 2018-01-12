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
library(tidyverse)

# Set up input directories.
pheno.dir = "/hpcdata/gac/derived/Attie_DO_Metabolomics/data/"
pheno.dict.dir = "/hpcdata/gac/raw/Attie_DO_Metabolomics/formatted_data/phenotype_dictionary/"
probs.dir = "/hpcdata/gac/derived/CGD_DO_Genoprobs/"
uwisc.dir = "/hpcdata/gac/raw/Attie_DO_Metabolomics/formatted_data/"
out.dir   = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/"

#########################################################################
# Cecum Metabolites: JAX normalization: sex, generation and batch covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_cecum_metabolites_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "cecum_metabolites_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$R_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_cecum_metabolites_zscore_normalized.rds"))
rownames(pheno.rz) = pheno.rz$Mouse.ID

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex    = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch  = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_cecum_metabolites_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#########################################################################
# Liver Metabolites: JAX normalization: sex, generation and batch covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "liver_metabolites_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_normalized.rds"))
rownames(pheno.rz) = pheno.rz$Mouse.ID

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex    = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch  = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_metabolites_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#################################################################
# Liver Metabolites: JAX normalization, sex and generation covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "liver_metabolites_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_normalized.rds"))
rownames(pheno.rz) = pheno.rz$Mouse.ID

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

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
pheno$Mouse.ID = sub("\\-|\\.", "", pheno$Mouse.ID)
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "liver_lipids_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_lipids_zscore_normalized.rds"))
pheno.rz$Mouse.ID = sub("\\-|\\.", "", pheno.rz$Mouse.ID)
rownames(pheno.rz) = pheno.rz$Mouse.ID
colnames(pheno.rz) = sub("Batch", "batch", colnames(pheno.rz))

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_liver_lipids_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


######################################################################
# Cecum Lipids: JAX normalization, sex, generation & batch covariates.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_cecum_lipids_normalized.rds"))
pheno$Mouse.ID = sub("\\-|\\.", "", pheno$Mouse.ID)
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "cecum_lipids_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_cecum_lipids_zscore_normalized.rds"))
pheno.rz$Mouse.ID = sub("\\-|\\.", "", pheno.rz$Mouse.ID)
rownames(pheno.rz) = pheno.rz$Mouse.ID

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch  = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_cecum_lipids_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)

#######################################################################
# Plasma Lipids: JAX normalization, sex, generation & batch covariates.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))
pheno$Mouse.ID = sub("\\-|\\.", "", pheno$Mouse.ID)
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "plasma_lipids_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_plasma_lipids_zscore_normalized.rds"))
pheno.rz$Mouse.ID = sub("\\-|\\.", "", pheno.rz$Mouse.ID)
rownames(pheno.rz) = pheno.rz$Mouse.ID

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_plasma_lipids_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


##########################################################################
# Metabolites: U. Wisc. normalization, sex, generation & batch covariates.

# Load in the phenotypes.
#pheno = read.delim(paste0(uwisc.dir, "14June2017_DOLiverMetabolites_NORM.txt"))
#rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
#covar = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))

#pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
#rownames(pheno) = pheno$Mouse.ID
#rm(covar)

# Create data dictionary.
#num.covar.columns = 13
#num.pheno.columns = ncol(pheno) - 13
#pheno.descr = data.frame(
#              name = colnames(pheno),
#              short_name = colnames(pheno),
#              description = colnames(pheno),
#              category = c(rep(NA, num.covar.columns), rep("metabolites", num.pheno.columns)),
#              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
#              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
#              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
#              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
#              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
#              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
#              )

# Load in Z-score phenotypes.
#pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
#genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
#samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
#pheno    = pheno[samples,]
#pheno.rz = pheno.rz[samples,]
#for(i in 1:length(genoprobs)) {
#  genoprobs[[i]] = genoprobs[[i]][samples,,]
#}

#print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
#pheno$sex = factor(pheno$sex)
#pheno$DOwave = factor(pheno$DOwave)
#pheno$batch = factor(pheno$batch)
#pheno.rz$sex = factor(pheno.rz$sex)
#pheno.rz$DOwave = factor(pheno.rz$DOwave)
#pheno.rz$batch = factor(pheno.rz$batch)

# Read in the marker data.
#markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
#map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
#K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
#save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
#     "attie_liver_metabolites_UWisc_qtl2_input.Rdata"))

#rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


##########################################################################
# Metabolites: U. Wisc. normalization, sex & generation covariates.

# Load in the phenotypes.
#pheno = read.delim(paste0(uwisc.dir, "14June2017_DOLiverMetabolites_NORM.txt"))
#rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
#covar = readRDS(paste0(pheno.dir, "attie_liver_metabolites_normalized.rds"))

#pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
#rownames(pheno) = pheno$Mouse.ID
#rm(covar)

# Create data dictionary.
#num.covar.columns = 13
#num.pheno.columns = ncol(pheno) - 13
#pheno.descr = data.frame(
#              name = colnames(pheno),
#              short_name = colnames(pheno),
#              description = colnames(pheno),
#              category = c(rep(NA, num.covar.columns), rep("metabolites", num.pheno.columns)),
#              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
#              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
#              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation", num.pheno.columns)),
#              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
#              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
#              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
#              )

# Load in Z-score phenotypes.
#pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_metabolites_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
#genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
#samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
#pheno    = pheno[samples,]
#pheno.rz = pheno.rz[samples,]
#for(i in 1:length(genoprobs)) {
#  genoprobs[[i]] = genoprobs[[i]][samples,,]
#}

#print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
#pheno$sex = factor(pheno$sex)
#pheno$DOwave = factor(pheno$DOwave)
#pheno$batch = factor(pheno$batch)
#pheno.rz$sex = factor(pheno.rz$sex)
#pheno.rz$DOwave = factor(pheno.rz$DOwave)
#pheno.rz$batch = factor(pheno.rz$batch)

# Read in the marker data.
#markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
#map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
#K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
#save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
#     "attie_liver_metabolites_UWisc_qtl2_input_sex_gen.Rdata"))

#rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#######################################
# Liver Lipids: U. Wisc. normalization.

# Load in the phenotypes.
#pheno = read.delim(paste0(uwisc.dir, "15June2017_DOLiverLipidomics.txt"))
#rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
#covar = readRDS(paste0(pheno.dir, "attie_liver_lipids_normalized.rds"))

#pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
#rownames(pheno) = pheno$Mouse.ID
#rm(covar)

# Create data dictionary.
#num.covar.columns = 13
#num.pheno.columns = ncol(pheno) - 13
#pheno.descr = data.frame(
#              name = colnames(pheno),
#              short_name = colnames(pheno),
#              description = colnames(pheno),
#              category = c(rep(NA, num.covar.columns), rep("lipids", num.pheno.columns)),
#              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
#              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
#              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
#              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
#              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
#              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
#              )

# Load in Z-score phenotypes.
#pheno.rz = readRDS(paste0(pheno.dir, "attie_liver_lipids_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
#genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
#samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
#pheno    = pheno[samples,]
#pheno.rz = pheno.rz[samples,]
#for(i in 1:length(genoprobs)) {
#  genoprobs[[i]] = genoprobs[[i]][samples,,]
#}

#print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
#pheno$sex = factor(pheno$sex)
#pheno$DOwave = factor(pheno$DOwave)
#pheno$batch = factor(pheno$batch)
#pheno.rz$sex = factor(pheno.rz$sex)
#pheno.rz$DOwave = factor(pheno.rz$DOwave)
#pheno.rz$batch = factor(pheno.rz$batch)

# Read in the marker data.
#markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
#map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
#K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
#save(pheno, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
#     "attie_liver_lipids_qtl2_UWisc_input.Rdata"))

#rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#########################################################################
# Plasma Metabolites: JAX normalization: sex, generation and batch covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_plasma_metabolites_normalized.rds"))
pheno$Mouse.ID = sub("\\-|\\.", "", pheno$Mouse.ID)
rownames(pheno) = pheno$Mouse.ID

pheno.descr = read.delim(paste0(pheno.dict.dir, "plasma_metabolites_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_plasma_metabolites_zscore_normalized.rds"))
pheno.rz$Mouse.ID = sub("\\-|\\.", "", pheno.rz$Mouse.ID)
rownames(pheno.rz) = pheno.rz$Mouse.ID

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex    = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch  = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
     "attie_plasma_metabolites_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)



#######################################
# Plasma Lipids: U. Wisc. normalization.

# Load in the phenotypes.
#pheno = read.delim(paste0(uwisc.dir, "15June2017_DOPlasmaLipidomics.txt"))
#rownames(pheno) = pheno$Mouse.ID

# Load in the covariates.
#covar = readRDS(paste0(pheno.dir, "attie_plasma_lipids_normalized.rds"))

#pheno = right_join(covar[,1:13], pheno, by = "Mouse.ID")
#rownames(pheno) = pheno$Mouse.ID
#rm(covar)

# Create data dictionary.
#num.covar.columns = 13
#num.pheno.columns = ncol(pheno) - 13
#pheno.descr = data.frame(
#              name = colnames(pheno),
#              short_name = colnames(pheno),
#              description = colnames(pheno),
#              category = c(rep(NA, num.covar.columns), rep("lipids", num.pheno.columns)),
#              is_numeric = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns)),
#              omit = c(rep(TRUE, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              units = c(rep(NA, num.covar.columns), rep("intensity", num.pheno.columns)),
#              is_derived = c(rep(NA, num.covar.columns), rep(FALSE, num.pheno.columns)),
#              covar_list = c(rep(NA, num.covar.columns), rep("sex:generation:Batch", num.pheno.columns)),
#              is_mouse_id = c(TRUE, rep(FALSE, ncol(pheno) - 1)),
#              is_covar = c(FALSE, rep(TRUE, 12), rep(TRUE, num.pheno.columns)),
#              is_pheno = c(rep(FALSE, num.covar.columns), rep(TRUE, num.pheno.columns))
#              )

# Load in Z-score phenotypes.
#pheno.rz = readRDS(paste0(pheno.dir, "attie_plasma_lipids_zscore_uwisc_normalized.rds"))

# Load in genoprobs.
#genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
#samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
#pheno    = pheno[samples,]
#pheno.rz = pheno.rz[samples,]
#for(i in 1:length(genoprobs)) {
#  genoprobs[[i]] = genoprobs[[i]][samples,,]
#}

#print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
#pheno$sex = factor(pheno$sex)
#pheno$DOwave = factor(pheno$DOwave)
#pheno$batch = factor(pheno$batch)
#pheno.rz$sex = factor(pheno.rz$sex)
#pheno.rz$DOwave = factor(pheno.rz$DOwave)
#pheno.rz$batch = factor(pheno.rz$batch)

# Read in the marker data.
#markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
#map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
#K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

# Save to *.Rdata file.
#save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, file = paste0(out.dir, 
#     "attie_plasma_lipids_qtl2_UWisc_input.Rdata"))

#rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)


#########################################################################
# Islet Proteins: JAX normalization: sex, generation and batch covar.

# Load in the phenotypes.
pheno = readRDS(paste0(pheno.dir, "attie_islet_proteins_normalized.rds"))
pheno$Mouse.ID = sub("\\-|\\.", "", pheno$Mouse.ID)
rownames(pheno) = pheno$Mouse.ID
colnames(pheno) = sub("^wave", "DOwave", colnames(pheno))
colnames(pheno) = sub("sex\\.x", "sex", colnames(pheno))
colnames(pheno) = sub("Batch", "batch", colnames(pheno))

pheno.descr = read.delim(paste0(pheno.dict.dir, "islet_proteins_pheno_dict.txt"))

stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(colnames(pheno) == pheno.descr$data_name)

# Load in Z-score phenotypes.
pheno.rz = readRDS(paste0(pheno.dir, "attie_islet_proteins_zscore_normalized.rds"))
pheno.rz$Mouse.ID = sub("\\-|\\.", "", pheno.rz$Mouse.ID)
rownames(pheno.rz) = pheno.rz$Mouse.ID
colnames(pheno.rz) = sub("^wave", "DOwave", colnames(pheno.rz))
colnames(pheno.rz) = sub("Batch", "batch", colnames(pheno.rz))

# Load in genoprobs.
genoprobs = readRDS("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/genoprobs/attie_DO500_genoprobs_v3.rds")

# Subset by samples.
samples  = sort(intersect(rownames(pheno), rownames(genoprobs[[1]])))
pheno    = pheno[samples,]
pheno.rz = pheno.rz[samples,]
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][samples,,]
}

print(paste(nrow(pheno), "Samples"))

# Create factors for covariates that we want to map as factors.
pheno$sex = factor(pheno$sex)
pheno$DOwave = factor(pheno$DOwave)
pheno$batch  = factor(pheno$batch)
pheno.rz$sex = factor(pheno.rz$sex)
pheno.rz$DOwave = factor(pheno.rz$DOwave)
pheno.rz$batch  = factor(pheno.rz$batch)

# Read in the marker data.
markers = readRDS("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_0.02cM_plus.rds")

# Make marker map.
map = map_df_to_list(map = markers[,1:4], pos_column = "pos")

# Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)

stopifnot(ncol(pheno) == ncol(pheno.rz))
stopifnot(ncol(pheno) == nrow(pheno.descr))
stopifnot(nrow(pheno) == nrow(genoprobs[[1]]))
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))

# Save to *.Rdata file.
save(pheno, pheno.rz, pheno.descr, genoprobs, K, map, markers, file = paste0(out.dir, 
     "attie_islet_proteins_qtl2_input.Rdata"))

rm(pheno, pheno.rz, pheno.descr, genoprobs, K, map)

