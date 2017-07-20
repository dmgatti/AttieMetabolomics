################################################################################
# Compare the four different metabolite mapping methods:
# JAX normalized: sex, gen & batch as covariates.
# U. Wisc. normalized: sex, gen & batch as covariates.
# JAX normalized: sex & gen as covariates.
# U. Wisc. normalized: sex & gen as covariates.
################################################################################
options(stringsAsFactors = F)
library(tidyverse)
library(broom)

setwd("/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/")
jax.dir = "metabolites_norm_jax/"
jax.sg.dir = "metabolites_norm_jax_sex_gen/"
uw.dir = "metabolites_norm_uwisc/"
uw.sg.dir = "metabolites_norm_uwisc_sex_gen/"

# Read in the QTL summaries with the maximum LOD for each analyte.
data = as.list(1:4)
names(data) = c("JAX norm: sex, gen & batch", "JAX norm: sex & gen", "U Wisc norm: sex, gen & batch", "U Wisc norm: sex & gen")
data[[1]] = read.csv(paste0(jax.dir, "liver_metabolites_jax_norm_qtl_summary.csv"))
data[[2]] = read.csv(paste0(jax.sg.dir, "liver_metabolites_jax_norm_sex_gen_qtl_summary.csv"))
data[[3]] = read.csv(paste0(uw.dir, "liver_metabolites_uwisc_norm_qtl_summary.csv"))
data[[4]] = read.csv(paste0(uw.sg.dir, "liver_metabolites_uwisc_norm_sex_gen_qtl_summary.csv"))

# Make empirical CDFs of the LOD scores for each data set.
plot_cdf = function(x, ...) {
  lines(sort(x), (1:length(x)) / length(x), ...)
}

max.lod = max(sapply(data, function(z) { max(z[,5]) }))

pdf("../../figures/JAX_Uwisc_metabolite_QTL_ecdfs.pdf")
plot(-1, -1, col = 0, xlim = c(4, 10), ylim = c(0, 1), xlab = "LOD", las = 1,
     ylab = "CDF", main = "Comparison of maximum LODs for metabolite using\ndifferent mapping and normalization methods")
abline(h = 0:5/5, col = "grey80")
abline(v = 4:10, col = "grey80")
for(i in 1:length(data)) {
  plot_cdf(x = data[[i]][,5], col = i, lwd = 2)
}
legend("bottomright", legend = names(data), lwd = 2, col = 1:4)
dev.off()


# Look at the correlation between analytes in the JAX and U Wisc data.
jax = readRDS("/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/data/attie_liver_metabolites_normalized.rds")
rownames(jax) = jax[,1]
jax = as.matrix(jax[,-(1:13)])
uw  = read.delim("/hpcdata/gac/raw/Attie_DO_Liver_Metabolomics/formatted_data/14June2017_DOLiverMetabolites_NORM.txt")
rownames(uw) = uw[,1]
uw = as.matrix(uw[,-1])

jax.cor = cor(jax)
uw.cor  = cor(uw)

pdf("../../figures/JAX_Uwisc_metabolite_cor.pdf")
plot(density(uw.cor), lwd = 2, las = 1, xlab = "Pearson Correlation", 
     main = "Correlation between Metabolites")
abline(h = 0:4/2,  col = "grey80")
abline(v = -2:2/2, col = "grey80")
lines(density(uw.cor),  lwd = 2, col = 1)
lines(density(jax.cor), lwd = 2, col = 2)
legend("topleft", col = 1:2, lwd = 2, legend = c("U. Wisc.", "JAX"), bg = "white")
dev.off()




