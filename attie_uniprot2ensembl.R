################################################################################
# Translate the UNIPROT IDs in the Attie Islet proteomics data to Ensembl
# protein IDs.
# Daniel Gatti
# dan.gatti@jax.org
# Sept. 27, 2017
################################################################################
options(stringsAsFacors = F)
library(tidyverse)
library(AnnotationHub)

uniprot.dir = "/hpcdata/gac/resource/Uniprot/"
input.dir = "/hpcdata/gac/derived/Attie_DO_Metabolomics/data/"

# Get Ensembl GTF
hub = AnnotationHub()
hub = query(hub, c("ensembl", "mus musculus", "gtf"))
ensembl = hub[[names(hub)[hub$title == "Mus_musculus.GRCm38.88.gtf"]]]
ensembl = ensembl[ensembl$type == "gene"]

# Read in the UNIPROT database.
uniprot = read_delim(file = paste0(uniprot.dir, "UP000000589_10090.idmapping.gz"), 
          delim = "\t", col_names = FALSE)

# Keep rows with gene ID, protein ID, Transcript ID and symbol.
dim(uniprot)
uniprot = uniprot %>% filter(X2 %in% c("Gene_Name", "STRING", "Ensembl", 
            "GeneTree"))
dim(uniprot)

# Spread the data up into one tibble per UNIPROT ID.
uniprot = split(uniprot, uniprot[[1]])

# Go through each entry and combine duplicates.
dupl = sapply(uniprot, function(z) { any(duplicated(z$X2)) })
dupl = which(dupl)

for(i in dupl) {

  wh = which(duplicated(uniprot[[i]]$X2))

} # for(i)



# Load in the Attie islet proteomics data.
prot = readRDS(paste0(input.dir, "attie_islet_proteins_normalized.rds"))

# Extract the column names and parse the UNIPROT IDs (_ delimited)
prot = colnames(prot)[-(1:13)]
prot = strsplit(prot, "_")

# Note: Gene_Name == symbol, Ensembl = ENSMUSG, GeneTree = ENSGT,
#       STRING == ENSMUSP
result = data.frame(uniprot = unlist(prot), Gene_Name = NA, Ensembl = NA, 
         GeneTree = NA, STRING = NA)

for(i in 1:length(prot)) {

  if(i %% 10 == 0) {
    print(paste(i, "of", length(prot)))
  } # if(i %% 10 == 0)

  m = match(prot[[i]], names(uniprot))
  m = m[!is.na(m)]

  if(length(m) > 0) {

    up = uniprot[m]

    for(j in 1:length(up)) {

      row = which(result$uniprot == names(up)[j])
      m2  = match(up[[j]]$X2, colnames(result))
      result[row, m2] = up[[j]]$X3

    } # for(j)

  } # if(any(!is.na(m)))

} # for(i)


