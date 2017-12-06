################################################################################
# For each association mapping file in the given directory, intersect them
# with the exons of protein coding genes.
#
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 29, 2017
################################################################################
options(stringsAsFactors = F)
library(rtracklayer)
library(AnnotationHub)
library(GenomicRanges)
library(qtl2)

# Arguments:
# arg1: assoc.dir: full path to the directory containing the association
#                  mapping files.
# arg2: lod.drop: LOD drop below the maximum.
# arg2
assoc.dir = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/"
lod.drop = 0.5
output.file = "/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_assoc_gene_intersection.txt"

args = commandArgs(trailingOnly = TRUE)
assoc.dir   = args[1]
lod.drop    = as.numeric(args[2])
output.file = args[3]

ensembl.file = "/hpcdata/gac/projects/Attie_DO_Metabolomics/Mus_musculus.GRCm38.90.gtf.rds"

##########################
# Get the Ensembl genes. #
##########################
ensembl = NULL
if(file.exists(ensembl.file)) {

  ensembl = readRDS(ensembl.file)

} else {

  hub = AnnotationHub()
  hub = AnnotationHub::query(hub, c("gtf", "mus musculus", "ensembl"))
  ensembl = hub[[names(hub)[hub$title == "Mus_musculus.GRCm38.90.gtf"]]]
  saveRDS(ensembl, file = ensembl.file)

} # else

ensembl = ensembl[ensembl$type == "exon"]

assoc.files = dir(path = assoc.dir, pattern = "_assoc_chr([0-9]+|X)",
              full.names = TRUE)

for(i in 1:length(assoc.files)) {

  print(paste(i, "of", length(assoc.files)))

  analyte = gsub(paste0("^", assoc.dir, "/|_assoc_chr([0-9]+|X).Rdata$"), "", assoc.files[i])

  # Load in 2 assoc objects: scan1output and snpinfo.
  load(assoc.files[i])
  
  map = qtl2plot:::snpinfo_to_map(assoc[[2]])
  tmp = qtl2plot:::expand_snp_results(assoc[[1]], map, assoc[[2]])
  assoc = data.frame(assoc[[2]], lod = tmp$lod[,1])
  rm(tmp, map)
  assoc = assoc[assoc$lod >= max(assoc$lod, na.rm = T) - lod.drop,]
  assoc = GRanges(seqnames = assoc$chr, ranges = IRanges(start = assoc$pos * 1e6,
          width = 1), alleles = assoc$alleles, sdp = assoc$sdp, csq = assoc$csq,
          lod = assoc$lod)

  # Intersect most signficant SNPs with genes.
  ol = findOverlaps(assoc, ensembl)

  if(length(ol) > 0) {

    hits = cbind(analyte, as.data.frame(assoc[queryHits(ol),]),
                 as.data.frame(ensembl[subjectHits(ol),]))
    colnames(hits) = make.unique(colnames(hits))

    hits = hits[,c("analyte", "seqnames", "start", "alleles", "sdp", "lod",
                "start.1", "end.1", "strand.1", "gene_id", "gene_name", 
                "gene_biotype", "transcript_id", "transcript_name",
                "exon_id", "protein_id")]

    write.table(hits, file = output.file, row.names = FALSE, quote = FALSE,
                sep = "\t", col.names = i == 1, append = i > 1)

  } # if(length(ol) > 0)

} # for(i)

