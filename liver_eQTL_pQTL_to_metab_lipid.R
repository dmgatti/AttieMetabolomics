################################################################################
# Compare the Svenson DO478 eQTL and pQTL peaks to the Attie DO liver 
# metabolomics.
#
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 21, 2017
################################################################################
options(stringsAsFactors = F)
library(GenomicRanges)
library(tidyverse)
library(qtl2)

source("/hpcdata/gac/projects/Attie_DO_Metabolomics/scripts/mqtl_eqtl_helper_fxns.R")

data.file = "/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_liver_metabolites_qtl2_input.Rdata"
summary.file = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_qtl_summary_thresh_6.csv"
mqtl.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/"
eqtl.dir = "/hpcdata/gac/derived/CGD_DO_Liver_RNASeq/QTL/Additive/"
out.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/eQTL_overlap/"
fig.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/"

thr = 7

# Load in the Svenson DO478 eQTL LOD table.
eqtl = read.csv("/projects/cgd/QTL_mapping/liver2/output/QTLAllAdditive.csv")
eqtl.gr = GRanges(seqnames = eqtl$chr, ranges = IRanges(start = eqtl$AdditivePos - 1e6, 
                  end = eqtl$AdditivePos + 1e6), ensembl = eqtl$id, 
                  symbol = eqtl$symbol, lod = eqtl$AdditiveLOD)
eqtl.gr = eqtl.gr[eqtl.gr$lod > thr]

# Load in the pQTL LOD table.


# Load in the data.
load(data.file)

# Load in the liver metabolite and lipid LOD table.
metab = read.csv("/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_qtl_summary_thresh_6.csv")
metab$chr = factor(metab$chr, levels = c(1:19, "X"))

png(paste0(fig.dir, "liver_metabolite_qtl_hotspots.png"), width = 1200, height = 800, 
    res = 128)
ggplot(metab, aes(pos)) +
  geom_histogram(binwidth = 5) +
  facet_grid(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(0.1, "lines")) + 
  labs(title = "Liver Metabolite QTL Hotspots")
dev.off()

metab.gr = GRanges(seqnames = metab$chr, ranges =
           IRanges(start = metab$pos * 1e6 - 1e6, end = metab$pos * 1e6 + 1e6),
           analyte = metab$analyte, lod = metab$lod)
metab.gr = metab.gr[metab.gr$lod > thr]

# Overlap QTL
metab.eqtl.ol = findOverlaps(metab.gr, eqtl.gr)

metab.eqtl.ol = cbind(as.data.frame(metab.gr[queryHits(metab.eqtl.ol)]), 
                      as.data.frame(eqtl.gr[subjectHits(metab.eqtl.ol)]))

metab.eqtl.ol = metab.eqtl.ol[,-grep("width|strand", colnames(metab.eqtl.ol))]

colnames(metab.eqtl.ol)[1] = "chr"

colnames(metab.eqtl.ol)[c(6:8, 11)] = c("gene_chr", "gene_start", "gene_end", "gene_lod")

metab.eqtl.ol$start = metab.eqtl.ol$start * 1e-6
metab.eqtl.ol$end   = metab.eqtl.ol$end   * 1e-6
metab.eqtl.ol$gene_start = metab.eqtl.ol$gene_start * 1e-6
metab.eqtl.ol$gene_end   = metab.eqtl.ol$gene_end * 1e-6

metab.eqtl.ol %>% group_by(analyte) %>%
  summarize(n = n()) %>%
  ggplot(aes(n)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Number of eQTL per Metabolite") 

spl = split(metab.eqtl.ol, metab.eqtl.ol$analyte)

result = NULL

for(i in 1:length(spl)) {

  ss = spl[[i]]

  for(j in 1:nrow(ss)) {

    filename = dir(path = eqtl.dir, pattern = paste0(ss$ensembl[j],
               "_coef_chr", ss$gene_chr[j], ".rds"), full.names = TRUE)

    if(length(filename) > 0) {

      # Get the metabolite coefficients and keep the ones within +/- 1 Mb of the peak.
      mqtl.coef = readRDS(paste0(mqtl.dir, ss$analyte[j], "_coef_chr", ss$chr[j], ".rds"))
      mqtl.coef = mqtl.coef[,LETTERS[1:8]]

      pos = as.numeric(sapply(strsplit(rownames(mqtl.coef), "_"), "[", 2)) * 1e-6

      mqtl.coef = mqtl.coef[pos >= ss$start[j] & pos <= ss$end[j],]
#  mqtl.coef = apply(mqtl.coef, 2, rank) / (nrow(mqtl.coef) + 1)
#  mqtl.coef = apply(mqtl.coef, 2, qnorm)

      gene.coef = readRDS(filename)[,LETTERS[1:8]]
      gene.coef = gene.coef[rownames(mqtl.coef),]
#      gene.coef = apply(gene.coef, 2, rank) / (nrow(gene.coef) + 1)
#      gene.coef = apply(gene.coef, 2, qnorm)
      cr = diag(cor(mqtl.coef, gene.coef))
      result = rbind(result, c(ss$analyte[j], ss$ensembl[j], ss$symbol[j], cr))

#layout(matrix(1:2, 2, 1))
#plot(mqtl.coef[,1], type = "l", col = CCcolors[1], ylim = range(mqtl.coef),
#     main = ss$analyte[j])
#for(k in 2:8) { lines(mqtl.coef[,k], col = CCcolors[k]) }
#plot(gene.coef[,1], type = "l", col = CCcolors[1], ylim = range(gene.coef),
#     main = ss$analyte[j])
#for(k in 2:8) { lines(gene.coef[,k], col = CCcolors[k]) }

    } else {

      warning(paste(ss$ensembl[i], "chr", ss$gene_chr[i], "coef file not found."))

    } # else

  } # for(j)

} # for(i)

result = data.frame(analyte = result[,1], 
                    ensembl = result[,2],
                    symbol  = result[,3],
                    A = as.numeric(result[,"A"]),
                    B = as.numeric(result[,"B"]),
                    C = as.numeric(result[,"C"]),
                    D = as.numeric(result[,"D"]),
                    E = as.numeric(result[,"E"]),
                    F = as.numeric(result[,"F"]),
                    G = as.numeric(result[,"G"]),
                    H = as.numeric(result[,"H"]))


