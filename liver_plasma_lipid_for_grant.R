# Mark requested plots for some lipids for a grant.
# In liver:
# analyte	marker	chr	pos	lod
# CE.18.1	JAX00546396	4	19.8	7.1
# CE.18.1	UNC17296803	9	122.8	6.9
# CE.18.1	UNC21853905	12	108.3	7.0
# CE.18.1	UNC24283981	14	76.6	6.2
# CE.22.6	UNCHS038832	14	77.2	6.4

# In plasma:
# analyte	marker	chr	pos	lod
# CE_20.5	JAX00012497	1	170.6	6.8
# CE_18.1	UNC2194776	1	171.1	9.1
# CE_18.2	UNC2194776	1	171.1	11.4
# CE_20.3	UNC2194776	1	171.1	9.2
# CE_20.4	UNC2194776	1	171.1	7.8
# CE_22.6	UNC2194776	1	171.1	8.8
# CE_20.5	JAX00569577	4	143.3	6.4
# CE_18.1	UNCHS031488	11	79.3	6.4
# CE_20.4	JAX00319703	11	99.3	6.2
# CE_20.5	UNCHS032029	11	100.1	6.1
# CE_20.3	UNC20262426	11	101.9	6.0
# CE_18.2	UNCHS047372	19	23.4	6.5
# CE_20.5	UNCHS047482	19	27.9	7.3

options(stringsAsFactors = FALSE)
library(qtl2)

liver.data.file  = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_liver_lipids_qtl2_input.Rdata"
plasma.data.file = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_plasma_lipids_qtl2_input.Rdata"
output.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/Mark_grant/"

liver = data.frame(analyte = c("CE.18.1", "CE.18.1", "CE.18.1", "CE.18.1", "CE.22.6"),
                   marker = c("JAX00546396", "UNC17296803", "UNC21853905", "UNC24283981", "UNCHS038832"),
                   chr = c(4, 9, 12, 14, 14),
                   pos = c(19.8, 122.8, 108.3, 76.6, 77.2),
                   lod = c(7.1, 6.9, 7, 6.2, 6.4))

plasma = data.frame(analyte = c("CE_20.5", "CE_18.1", "CE_18.2", "CE_20.3",
                    "CE_20.4", "CE_22.6", "CE_20.5", "CE_18.1", "CE_20.4",
                    "CE_20.5", "CE_20.3", "CE_18.2", "CE_20.5"),
                    marker = c("JAX00012497", "UNC2194776", "UNC2194776", "UNC2194776", "UNC2194776",
                    "UNC2194776", "JAX00569577", "UNCHS031488", 
                    "JAX00319703", "UNCHS032029", "UNC20262426", 
                    "UNCHS047372", "UNCHS047482"),
                    chr = c(1, 1, 1, 1, 1, 1, 4, 11, 11, 11, 11, 19, 19),
                    pos = c(170.6, 171.1, 171.1, 171.1, 171.1, 171.1, 143.3, 79.3, 99.3, 100.1, 101.9, 23.4, 27.9),
                    lod = c(6.8, 9.1, 11.4, 9.2, 7.8, 8.8, 6.4, 6.4, 6.2, 6.1, 6, 6.5, 7.3))

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

#########
# Liver #
#########
# Read in the liver data.
load(liver.data.file)

# Make covariates.
covar = model.matrix(~sex + DOwave + batch, data = pheno)[,-1]

for(i in 1:nrow(liver)) {

  analyte = liver$analyte[i]
  chr     = liver$chr[i]

  index = which(colnames(pheno) == analyte)
  ph = pheno[,index,drop = FALSE]
#  ph[,1] = rankZ(ph[,1])

  qtl = scan1(genoprobs = genoprobs, pheno = ph, kinship = K,
              addcovar = covar, cores = 5)

  saveRDS(qtl, file = paste0(output.dir, "liver_", analyte, "_QTL.rds"))

  pdf(paste0(output.dir, "liver_", analyte, "_QTL.pdf"), width = 8, height = 6)
  plot_scan1(qtl, map, main = analyte)
  dev.off()

  blup = scan1blup(genoprobs = genoprobs[,chr], pheno = ph, 
                   kinship = K[chr], addcovar = covar, cores = 5,
                   quiet = FALSE)

  saveRDS(blup, file = paste0(output.dir, "liver_", analyte, "_coef_chr",
          chr, ".rds"))

  pdf(paste0(output.dir, "liver_", analyte, "_coef_chr", chr, ".pdf"), 
      width = 8, height = 6)
  plot_coefCC(blup, map, main = analyte, scan1_output = qtl)
  dev.off()

} # for(i)

rm(pheno, pheno.rz, covar, genoprobs, K)


##########
# Plasma #
##########
# Read in the liver data.
load(plasma.data.file)

# Make covariates.
covar = model.matrix(~sex + DOwave + batch, data = pheno)[,-1]

for(i in 1:nrow(plasma)) {

  analyte = plasma$analyte[i]
  chr     = plasma$chr[i]

  index = which(colnames(pheno) == analyte)
  ph = pheno[,index,drop = FALSE]

  qtl = scan1(genoprobs = genoprobs, pheno = ph, kinship = K,
              addcovar = covar, cores = 5)

  saveRDS(qtl, file = paste0(output.dir, "plasma_", analyte, "_QTL.rds"))

  pdf(paste0(output.dir, "plasma_", analyte, "_QTL.pdf"), width = 8, height = 6)
  plot_scan1(qtl, map, main = analyte)
  dev.off()

  blup = scan1blup(genoprobs = genoprobs[,chr], pheno = ph, 
                   kinship = K[chr], addcovar = covar, cores = 5,
                   quiet = FALSE)

  saveRDS(blup, file = paste0(output.dir, "plasma_", analyte, "_coef_chr",
          chr, ".rds"))

  pdf(paste0(output.dir, "plasma_", analyte, "_coef_chr", chr, ".pdf"), 
      width = 8, height = 6)
  plot_coefCC(blup, map, main = analyte, scan1_output = qtl)
  dev.off()

} # for(i)


