library(tidyverse)
source("/projects/dgatti/scripts/transcriptome_map.R")

data = read_csv("/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl/pQTL_AnalyteFamilies_UniProt_mappedToMGI_geneSymbol_geneLocation_withCisInfo.csv")

# Change the column names.
data = data %>% rename_all(tolower)
data = data %>% rename(lod = lod_qtl)
data$start_gene = data$start_gene * 1e-6
data$end_gene = data$end_gene * 1e-6
data$cis = data$cis == "True"

pdf("islet_protein_tmap.pdf", width = 8.5, height = 11)
tmap(data)
dev.off()
