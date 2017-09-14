################################################################################
# Normalize and impute missing data in the Attie islet protein data set.
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

# Read in the raw metabolomics data.
prot = read_delim(paste0(input.dir, "formatted_data/DO_islet_proteomics_non_normalized.txt"),
        delim = "\t")
rownames(prot) = prot$Mouse.ID

# Read in the sample annotation.
annot = read_delim(paste0(input.dir, "attie_DO_sample_annot.txt"), delim = "\t")

# Merge the sample annotation and data.
prot = right_join(annot, prot, by = "Mouse.ID")

# Keep columns that are < 50% NA.
dim(prot)
prot = prot[,colSums(is.na(prot)) < 0.5 * nrow(prot)]
dim(prot)

# Split up the sample annotation from the data and convert the data into a 
# numeric matrix.
annot = as.data.frame(prot[,1:10])
data  = as.matrix(prot[,-(1:10)])
rownames(data)  = annot$Mouse.ID

dim(data)

# 439 samples, 5433 proteins.

# Make a PCA plot of all of the data, with sample labels.
pc.data = pca(log(data), method = "bpca", nPcs = 5)
pdf("figures/islet_proteins_unnormalized_all_data_PCA.pdf")
batch.colors = as.numeric(factor(annot$Injection_batch))
plot(scores(pc.data), pch = 16, col = 0, main = "Un-normalized Metabolites, Colored by Batch")
text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data), 
     col = batch.colors)
dev.off()

# Remove control samples.
ctrl = grep("Std", annot$Mouse.ID)
data = data[-ctrl,]
annot = annot[-ctrl,]

# 375 samples (still with some replicates) and 5,433 analytes.
dim(data)
length(unique(rownames(data)))
rownames(data)[duplicated(rownames(data))]

######################
# Impute missing data.
data.log = log(data)

# pcaMethods wants samples in rows and variables in columns.
pc.data = pca(data.log, method = "bpca", nPcs = 5)
plot(pc.data)
abline(h = 0.95, col = 2)

# Make PCA plots of the unnormalized data, colored by batch, sex, etc.
pdf("figures/islet_proteins_unnormalized_PCA.pdf")

sex = factor(annot$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Un-normalized Metabolites Colored by Sex")
legend("topright", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

batch = factor(annot$Injection_batch)
batch.colors = rainbow(length(levels(batch)) - 1)[batch]
plot(scores(pc.data), pch = 16, col = batch.colors,
     main = "Un-normalized Metabolites Colored by Batch")
legend("topright", legend = levels(batch), pch = 16, col = batch.colors,
       y.intersp = 0.7)

wave = factor(annot$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Un-normalized Metabolites Colored by Wave")
legend("topright", legend = levels(wave), pch = 16, col = 1:length(levels(wave)))

diet.days = factor(annot$diet_days, levels = sort(unique(annot$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Un-normalized Metabolites Colored by Diet Days")

plate = factor(annot$Plate_number, levels = sort(unique(annot$Plate_number)))
plate.colors = rainbow(length(levels(plate)) - 1)[plate]
plot(scores(pc.data), pch = 16, col = plate.colors,
     main = "Un-normalized Metabolites Colored by Plate")
legend("topright", legend = levels(plate), pch = 16, col = plate.colors)

dev.off()

# Set up batch and model for comBat.
mod = model.matrix(~sex, data = annot)[,-3]
batch = annot$Injection_batch

chg = 1e6
iter = 1
repeat( {

  print(paste("Iteration", iter))

  # Impute missing data.
  miss = which(is.na(data.log))
  print(paste(length(miss), "missing points."))
  pc.data = pca(data.log, method = "bpca", nPcs = 7)
  data.compl = completeObs(pc.data)

  # Batch adjust.
  # ComBat wants the data with variable in rows and samples in columns.
  data.cb = ComBat(dat = t(data.compl), batch = batch, mod = mod)
  data.cb = t(data.cb)

  # Calculate the change.
  chg = sum((data.compl[miss] - data.cb[miss])^2)
  print(paste("   SS Change:", chg))

  # Put the missing data back in an impute again.
  if(chg > 10 & iter < 20) {

    data.cb[miss] = NA
    data.log = data.cb
    iter = iter + 1    

  } else {

    data.log = data.cb
    break
  }}
)

# Remove duplicate samples, retaining the one with the least missing data.
dupl = which(duplicated(rownames(data.log)))
dupl.data = data.log[rownames(data.log) %in% rownames(data.log)[dupl],]

# Keep the sample with the lowest no-call rate.
stopifnot(rownames(data) == rownames(data.log))

prop.missing = rowMeans(is.na(data))

unique.samples = unique(rownames(data.log))
keep = rep(FALSE, nrow(data.log))
for(i in 1:length(unique.samples)) {

  sample = unique.samples[i]
  wh = which(rownames(data.log) == sample)
  wh = wh[which.min(prop.missing[wh])]
  keep[wh] = TRUE

} # for(i)

data.log = data.log[keep,]
annot = annot[match(rownames(data.log), annot$Mouse.ID),]

# Merge in the Chr M and Y info.
attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
annot = right_join(annot, attie_MY, by = "Mouse.ID")
annot = annot[,c(1:10, 13:15)]
colnames(annot) = sub("\\.x", "", colnames(annot))

data.log = data.frame(Mouse.ID = rownames(data.log), data.log)
data.out = right_join(annot, data.log, by = "Mouse.ID")

saveRDS(data.out, file = paste0(output.dir, "attie_islet_proteins_normalized.rds"))

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

for(i in 14:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}

saveRDS(data.rz, file = paste0(output.dir, "attie_islet_proteins_zscore_normalized.rds"))


# Make PCA plots of the normalized data, colored by batch, sex, etc.
pdf("figures/metabolites_normalized_PCA.pdf", width = 12, height = 7)

pc.data = pca(as.matrix(data.out[,-(1:13)]), method = "bpca", nPcs = 5)

layout(matrix(1:2, 1, 2))
sex = factor(data.out$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Normalized Metabolites Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "Normalized Metabolites Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(data.out$Injection_batch)
batch.colors = rainbow(length(levels(batch)))[batch]
plot(scores(pc.data), pch = 16, col = batch.colors,
     main = "Normalized Metabolites Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = batch.colors),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = batch.colors,
     main = "Normalized Metabolites Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = batch.colors),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(data.out$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Normalized Metabolites Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "Normalized Metabolites Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Normalized Metabolites Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "Normalized Metabolites Colored by Diet Days")

plate = factor(annot$Plate_number, levels = sort(unique(annot$Plate_number)))
plate.colors = rainbow(length(levels(plate)) - 1)[plate]
plot(scores(pc.data), pch = 16, col = plate.colors,
     main = "Normalized Metabolites Colored by Plate")
legend("bottomleft", legend = levels(wave), pch = 16, col = plate.colors,
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = plate.colors,
     main = "Normalized Metabolites Colored by Plate")
legend("bottomleft", legend = levels(wave), pch = 16, col = plate.colors,
       x.intersp = 0.7, y.intersp = 0.7)

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
annot = data.out[,1:13]
data  = as.matrix(data.out[,-(1:13)])

pdf("figures/metabolites_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(data, range = 0)
dev.off()

pdf("figures/metabolites_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot$Injection_batch))]
heatmap(data, RowSideColors = batch.colors)
dev.off()



###########################
# Compare to U. Wisc. data.
norm = read.delim("/hpcdata/gac/raw/Attie_DO_Liver_Metabolomics/formatted_data/14June2017_DOLiverMetabolites_NORM.txt")
rownames(norm) = norm$Mouse.ID

# Transform each analyte into Z-scores.
norm.rz = norm

for(i in 2:ncol(norm.rz)) {
  norm.rz[,i] = rankZ(norm.rz[,i])
}

saveRDS(norm.rz, file = paste0(output.dir, "attie_islet_proteins_zscore_uwisc_normalized.rds"))

norm = as.matrix(norm[,-1])

# Merge the sample annotation with the U. Wisc. normalized data.
annot.wisc = annot[annot$Mouse.ID %in% rownames(norm),]
norm = norm[annot.wisc$Mouse.ID,]
stopifnot(annot.wisc$Mouse.ID == rownames(norm))

pc.data = pca(norm, method = "bpca", nPcs = 5)

pdf("figures/metabolites_UWisc_normalized_PCA.pdf", width = 12, height = 7)

layout(matrix(1:2, 1, 2))
sex = factor(annot.wisc$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Metabolites Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Metabolites Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(annot.wisc$Injection_batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Metabolites Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Metabolites Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(annot.wisc$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Metabolites Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Metabolites Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(annot.wisc$diet_days, levels = sort(unique(annot.wisc$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Metabolites Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Metabolites Colored by Diet Days")

dev.off()


# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
pdf("figures/metabolites_UWisc_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(norm, range = 0)
dev.off()

pdf("figures/metabolites_UWisc_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot.wisc$Injection_batch))]
heatmap(norm, RowSideColors = batch.colors)
dev.off()

