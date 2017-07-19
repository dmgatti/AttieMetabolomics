# Quick look at metabolite data.
options(stringsAsFactors = F)
library(tidyverse)
library(sva)
library(pcaMethods)

data.dir = "/hpcdata/gac/raw/Attie_DO_Liver_Metabolomics/formatted_data/"

working.dir = "/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/"

setwd(working.dir)

# Read in Attie sample meta-data.
covar = read.csv("/hpcdata/gac/raw/Attie_DO_Liver_Metabolomics/attie_sample_info.csv")



# Raw metabolites.
met = read_delim(paste0(data.dir, "7July2017_DOLiverMetabolites_RAW.txt"),
      delim = "\t")
colnames(met)[-(1:5)] = paste0("RT", colnames(met)[-(1:5)])
do.samples = grep("^[0-9]+$", met$DO.Number)
met$DO.Number[do.samples] = paste0("DO-", met$DO.Number[do.samples])
colnames(met)[3] = "Mouse.ID"

# Fix sample IDs.
wh = which(nchar(met$Mouse.ID) == 5)
met$Mouse.ID[wh] = sub("^DO-", "DO-0", met$Mouse.ID[wh])

# Image of the pattern of missing data.
tmp = as.matrix(met[,-(1:5)])
tmp = log(tmp, 2)
rownames(tmp) = met$Mouse.ID

batch = as.numeric(factor(met$Batch))

png("figures/metabolites_missing_data.png", width = 1000, height = 2000, res = 128)
layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.99), mgp = c(2, 0.5, 0), las = 1)
image(1:nrow(tmp), 1:ncol(tmp), tmp == 0, breaks = c(-0.5, 0.5, 1.5),
      col = c("white", "black"), xlab = "Samples", ylab = "Analytes")
par(plt = c(0.1, 0.99, 0.2, 0.99), mgp = c(2, 0.5, 0))
plot(rowMeans(tmp == 0), pch = 16, col = batch, xlab = "Samples", ylab = "Proportion == 0",
     xaxs = "i")
dev.off()

# Plot of mean intensity versus prop. missing.
png("figures/metabolite_mean_signal_vs_missing.png", width = 800, height = 800,
   res = 128)
mean.sig = colMeans(tmp, na.rm = T)
p.miss = colMeans(tmp == 0) * 100
plot(mean.sig, p.miss, pch = 16, col = rgb(0,0,0,0.3), xlab = "Mean log2(signal)",
     ylab = "% missing data", las = 1)
mod = lm(p.miss ~ mean.sig)
abline(mod, col = 2)
dev.off()

# PCA plot of raw data.
met.num = as.matrix(met[,-(1:7)])
met.num[met.num == 1] = NA
met.num = log(met.num)
met.scl = apply(met.num, 2, scale)
pc = pca(met.scl, nPcs = 5, method = "ppca")

png("figures/metabolites_PCA.png", width = 800, height = 800, res = 128)
plot(scores(pc), pch = 16, col = as.numeric(factor(met$Batch)), las = 1)
dev.off()

# There are 2 bad samples: DO-109 and a control sample.
# Remove them (rows 44 & 151)
remove = which(scores(pc)[,1] > 80)

met = met[-remove,]
met.num = met.num[-remove,]
met.scl = met.scl[-remove,]
pc = pca(met.scl, nPcs = 3, method = "ppca")
plot(scores(pc), pch = 16, col = as.numeric(factor(met$Batch)))

# Plot data by sample ID and batch
layout(matrix(1:2, 1, 2))
plot(scores(pc), col = 0, main = "Labeled by Sample ID")
text(scores(pc)[,1], scores(pc)[,2], labels = met$Mouse.ID)
plot(scores(pc), col = 0, main = "Labeled by Batch")
text(scores(pc)[,1], scores(pc)[,2], labels = sub("^Batch", "", met$Batch),
     col = as.numeric(factor(met$Batch)))

# The controls aren't helping us. Remove them.
met = met[met$Mouse.ID != "Control",]

# Merge in the known covariates.
colnames(covar)[1] = "Mouse.ID"
met = right_join(covar, met, by = "Mouse.ID")

# Normalize by batch. Use only batch 1.
batch = factor(met$Batch)
mod = model.matrix(~sex, data = met)

met.num = as.matrix(met[,-(1:7)])
met.num = log(met.num)

met.cb = ComBat(dat = t(met.num), batch = batch, mod = mod, par.prior = TRUE, 
             prior.plots = TRUE)
met.cb = cbind(met[,1:7], t(met.cb))

pc = pca(met.cb, nPcs = 3, method = "ppca")
plot(scores(pc), pch = sub("^Batch", "", met$Batch), col = as.numeric(factor(met$Batch)))
plot(scores(pc), pch = 16, col = as.numeric(factor(met$Que.Number)))
plot(scores(pc), pch = 16, col = as.numeric(factor(met$gen)))

plot(scores(pc), pch = 16, col = as.numeric(factor(met$sex)),
     main = "Colored by Sex")

pdf("figures/met_boxplot_by_batch.pdf")
for(i in 2:309) {
  print(met %>% select(matches("^(RT|Bat)")) %>%
    ggplot(aes_string(x = "Batch", y = colnames(.)[i])) + 
    geom_boxplot() + 
    scale_y_log10() +
    labs(title = colnames(met)[i+4]))
}
dev.off()

pdf("figures/met_boxplot_by_batch_after_combat.pdf")
for(i in 2:309) {
  print(met.cb %>% select(matches("^(RT|Bat)")) %>%
    ggplot(aes_string(x = "Batch", y = colnames(.)[i])) + 
    geom_boxplot() + 
    scale_y_log10() +
    labs(title = colnames(met.cb)[i+6]))
}
dev.off()

# Try setting the '1' values to NA and imputing.
met.num[met.num < 1e-8] = NA

pc = pca(met.num, nPcs = 3, method = "ppca")
met.imp = completeObs(pc)

met.imp.cb = ComBat(dat = t(met.imp), batch = batch, mod = mod, par.prior = TRUE, 
             prior.plots = TRUE)
met.imp.cb = cbind(met[,1:7], t(met.imp.cb))

pc = pca(met.imp.cb, nPcs = 3, method = "ppca")
plot(scores(pc), pch = 16, col = as.numeric(factor(met$sex)),
     main = "Colored by Sex")

plot(scores(pc), pch = 16, col = as.numeric(factor(met$Batch)),
     main = "Colored by Batch")


# Normalized metabolites from Univ. Wisc.
met2 = read_delim(paste0(data.dir, "14June2017_DOLive Metabolites_NORM.txt"),
       delim = "\t")
met2$Mouse.ID = paste0("DO-", met2$Mouse.ID)
wh = which(nchar(met2$Mouse.ID) == 5)
met2$Mouse.ID[wh] = sub("^DO-", "DO-0", met2$Mouse.ID[wh])
met2 = right_join(covar, met2, by = "Mouse.ID")

met2.num = as.matrix(met2[,-(1:3)])
met2.num.scl = apply(met2.num, 2, scale)

pc = pca(met2.num.scl, nPcs = 3, method = "ppca")
plot(scores(pc), pch = 16, col = as.numeric(factor(met2$sex)),
     main = "Colored by Sex")

plot(scores(pc), pch = 16, col = as.numeric(factor(met2$gen)),
     main = "Colored by Generation")


# Lipid data.
# Raw data.
lip1 = read_delim(paste0(data.dir, "21June2017_DOLiverLipidomicsRawMPK.txt"),
       delim = "\t")
wh = which(nchar(lip1$Mouse.ID) == 2)
lip1$Mouse.ID[wh] = paste0("0", lip1$Mouse.ID[wh])
lip1$Mouse.ID = paste0("DO-", lip1$Mouse.ID)
lip1$Mouse.ID = sub("DO-Control", "Control", lip1$Mouse.ID)

# Normalized data (by U. Wisc.)
lip2 = read_delim(paste0(data.dir, "15June2017_DOLiverLipidomics.txt"),
       delim = "\t")

# Image of the pattern of missing data.
tmp = as.matrix(lip1[,-(1:5)])
tmp = log(tmp + 1, 2)
rownames(tmp) = lip1$Mouse.ID

batch = as.numeric(factor(met$Batch))

png("figures/lipids_missing_data.png", width = 1000, height = 2000, res = 128)
layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.99), mgp = c(2, 0.5, 0), las = 1)
image(1:nrow(tmp), 1:ncol(tmp), tmp < 1, breaks = c(-0.5, 0.5, 1.5),
      col = c("white", "black"), xlab = "Samples", ylab = "Analytes")
par(plt = c(0.1, 0.99, 0.2, 0.99), mgp = c(2, 0.5, 0))
plot(rowMeans(tmp == 0), pch = 16, col = batch, xlab = "Samples", ylab = "Proportion == 0",
     xaxs = "i")
dev.off()

# Plot of mean intensity versus prop. missing.
png("figures/metabolite_mean_signal_vs_missing.png", width = 800, height = 800,
   res = 128)
mean.sig = colMeans(tmp, na.rm = T)
p.miss = colMeans(tmp < 1) * 100
plot(mean.sig, p.miss, pch = 16, col = rgb(0,0,0,0.3), xlab = "Mean log2(signal)",
     ylab = "% missing data", las = 1)
mod = lm(p.miss ~ mean.sig)
abline(mod, col = 2)
dev.off()

# PCA plot of raw data.
lip.num = as.matrix(lip1[,-(1:5)])
lip.num[lip.num == 0] = NA
lip.num = log(lip.num)
lip.scl = apply(lip.num, 2, scale)
pc = pca(lip.scl, nPcs = 5, method = "ppca")

png("figures/lipids_PCA.png", width = 800, height = 800, res = 128)
plot(scores(pc), pch = 16, col = 0, las = 1)
text(scores(pc)[,1], scores(pc)[,2], labels = lip1$Mouse.ID, 
     col = as.numeric(factor(lip1$Batch)))
dev.off()

plot(scores(pc), pch = 16, col = 0, las = 1)
text(scores(pc)[,1], scores(pc)[,2], labels = lip1$Mouse.ID, 
     col = as.numeric(factor(lip1$sex)))

lip1 = lip1[lip1$Mouse.ID != "Control",]
colnames(covar)[1] = "Mouse.ID"
lip1 = right_join(covar, lip1, by = "Mouse.ID")

# Combat normalize.
batch = factor(lip1$Batch)
mod = model.matrix(~sex, data = lip1)

lip1.num = as.matrix(lip1[,-(1:9)])
lip1.num = log(lip1.num)

lip1.cb = ComBat(dat = t(lip1.num), batch = batch, mod = mod, par.prior = TRUE, 
                 prior.plots = TRUE)
lip1.cb = cbind(lip1[,1:9], t(lip1.cb))

lip1.num = t(lip1.cb)
lip1.num.scl = apply(lip1.num, 2, scale)

pc = pca(lip1.num.scl, nPcs = 5, method = "ppca")

png("figures/lipids_PCA_after_batch_norm.png", width = 800, height = 800, 
    res = 128)
plot(scores(pc), pch = 16, col = 0, las = 1)
text(scores(pc)[,1], scores(pc)[,2], labels = lip1$Mouse.ID, 
     col = as.numeric(factor(lip1$Batch)))
dev.off()

png("figures/lipids_PCA_after_batch_norm_by_sex.png", width = 800, height = 800, 
    res = 128)
plot(scores(pc), pch = 16, col = 0, las = 1)
text(scores(pc)[,1], scores(pc)[,2], labels = lip1$Mouse.ID, 
     col = as.numeric(factor(lip1$sex)))
dev.off()
