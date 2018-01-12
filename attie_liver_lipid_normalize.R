################################################################################
# Normalize and impute missing data in the Attie liver lipid data set.
# Daniel Gatti
# dan.gatti@jax.org
# July 14, 2017
################################################################################
options(stringsAsFactors = F)
library(tidyverse)
library(pcaMethods)
library(sva)

input.dir  = "/hpcdata/gac/raw/Attie_DO_Metabolomics/"
output.dir = "/hpcdata/gac/derived/Attie_DO_Metabolomics/data/"

setwd("/hpcdata/gac/projects/Attie_DO_Metabolomics/")

# Read in the raw lipid data.
lipid = read_delim(paste0(input.dir, "formatted_data/03_January_2018_DO_Liver_Lipidomics_Raw.txt"),
        delim = "\t")

# Read in the sample annotation.
annot = read_delim(paste0(input.dir, "attie_DO_sample_annot.txt"), delim = "\t")
annot$Mouse.ID = gsub("[^[:alnum:]]", "", annot$Mouse.ID)

# Merge the sample annotation and data.
lipid = right_join(annot, lipid, by = "Mouse.ID")

# Split up the sample annotation from the data and convert the data into a 
# numeric matrix.
annot = as.data.frame(lipid[,1:8])
data  = as.matrix(lipid[,-(1:8)])
rownames(data)  = annot$Mouse.ID

colnames(annot) = sub("Batch", "batch", colnames(annot))

dim(data)
range(data)

# Make a PCA plot of all of the data, with sample labels.
pc.data = pca(log(data), method = "bpca", nPcs = 10)

pdf("figures/liver_lipids_unnormalized_all_data_PCA.pdf")

batch.colors = as.numeric(factor(annot$batch))
plot(scores(pc.data), pch = 16, col = 0, main = "Un-normalized Liver Lipids, Colored by Batch")
text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data), 
     col = batch.colors)

dev.off()

# Remove control samples.
ctrl  = grep("Control", annot$Mouse.ID)
data  = data[-ctrl,]
annot = annot[-ctrl,]

rownames(annot) = annot$Mouse.ID

stopifnot(rownames(annot) == rownames(data))

# 384 samples and 1562 analytes.
dim(data)

######################
# Impute missing data.
data.log = log(data)

# pcaMethods wants samples in rows and variables in columns.
pc.data = pca(data.log, method = "bpca", nPcs = 10)
plot(pc.data)
abline(h = 0.95, col = 2)

# Make PCA plots of the unnormalized data, colored by batch, sex, etc.
pdf("figures/liver_lipids_unnormalized_PCA.pdf")

sex = factor(annot$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Un-normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

batch = factor(annot$batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Un-normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       y.intersp = 0.7)

wave = factor(annot$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Un-normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)))

diet.days = factor(annot$diet_days, levels = sort(unique(annot$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Un-normalized Liver Lipids Colored by Diet Days")

dev.off()

boxplot(scores(pc.data)[,1] ~ batch)

# Set up batch and model for comBat.
annot$sex  = factor(annot$sex)
annot$wave = factor(annot$wave)
mod = model.matrix(~sex, data = annot)
batch = annot$batch

# Batch adjust.
# ComBat wants the data with variable in rows and samples in columns.
data.cb = ComBat(dat = t(data.log), batch = batch, mod = mod, prior.plots = TRUE)
data.cb = t(data.cb)

# No duplicate samples.
dupl = which(duplicated(rownames(data.cb)))

# Merge in the Chr M and Y info.
attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
attie_MY$Mouse.ID = gsub("[^[:alnum:]]", "", attie_MY$Mouse.ID)
annot = right_join(annot, attie_MY, by = "Mouse.ID")
annot = annot[,-grep("\\.y", colnames(annot))]
colnames(annot) = sub("\\.x", "", colnames(annot))

data.cb  = data.frame(Mouse.ID = rownames(data.cb), data.cb)
data.out = right_join(annot, data.cb, by = "Mouse.ID")
rownames(data.out) = data.out$Mouse.ID
colnames(data.out) = sub("wave", "DOwave", colnames(data.out))

saveRDS(data.out, file = paste0(output.dir, "attie_liver_lipids_normalized.rds"))

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

for(i in 11:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}

saveRDS(data.rz, file = paste0(output.dir, "attie_liver_lipids_zscore_normalized.rds"))


# Make PCA plots of the normalized data, colored by batch, sex, etc.
pdf("figures/liver_lipids_normalized_PCA.pdf", width = 12, height = 7)

pc.data = pca(as.matrix(data.out[,-(1:13)]), method = "bpca", nPcs = 10)

layout(matrix(1:2, 1, 2))
sex = factor(data.out$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(data.out$batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(data.out$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Normalized Liver Lipids Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "Normalized Liver Lipids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
annot = data.out[,1:11]
data  = as.matrix(data.out[,-(1:11)])

pdf("figures/liver_lipids_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(data, range = 0)
dev.off()

pdf("figures/liver_lipids_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot$batch))]
heatmap(data, RowSideColors = batch.colors)
dev.off()


###########################
# Compare to U. Wisc. data.
norm = read.delim("/hpcdata/gac/raw/Attie_DO_Metabolomics/formatted_data/15June2017_DOLiverLipidomics.txt")
rownames(norm) = norm$Mouse.ID

# Transform each analyte into Z-scores.
norm.rz = norm

for(i in 2:ncol(norm.rz)) {
  norm.rz[,i] = rankZ(norm.rz[,i])
}

saveRDS(norm.rz, file = paste0(output.dir, "attie_liver_lipids_zscore_uwisc_normalized.rds"))

norm = as.matrix(norm[,-1])

# Merge the sample annotation with the U. Wisc. normalized data.
annot.wisc = annot[annot$Mouse.ID %in% rownames(norm),]
norm = norm[annot.wisc$Mouse.ID,]
stopifnot(annot.wisc$Mouse.ID == rownames(norm))

pc.data = pca(norm, method = "bpca", nPcs = 10)

pdf("figures/liver_lipids_UWisc_normalized_PCA.pdf", width = 12, height = 7)

layout(matrix(1:2, 1, 2))
sex = factor(annot.wisc$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(annot.wisc$batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(annot.wisc$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(annot.wisc$diet_days, levels = sort(unique(annot.wisc$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Liver Lipids Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Liver Lipids Colored by Diet Days")

dev.off()


# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
pdf("figures/liver_lipids_UWisc_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(norm, range = 0)
dev.off()

pdf("figures/liver_lipids_UWisc_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot.wisc$batch))]
heatmap(norm, RowSideColors = batch.colors)
dev.off()

