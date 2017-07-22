########################
# Association mapping. #
########################
# NOTE: You need /data/gac/resource/CCsnps/ccfoundersnps.sqlite, which is 
#        a database of Sanger SNPs for the 8 DO founders.

# Arguments:
# probs: genoprobs object in qtl2 format.
# pheno: matrix of phenotypes, samples in rows, phenotypes in columns.
# addcovar: covariates matrix as used in scan1().
# intcovar: covariate to interact with QTL.
# k: list of kinship matrices.
# chr: Chromosome to map on.
# start: start position for mapping in Mb.
# end: end position for mapping in Mb.
# ncores: number of cores to use in mappin.
# map: list of map positions from qtl2.
# db.file: Location of the mySQL database containing the Sanger SNPs.
assoc_map = function(chr, start, end, probs, pheno, K, addcovar, intcovar = NULL, map,
            db.file = "/data/gac/resource/CCsnps/ccfoundersnps.sqlite", ncores) {

  my_db = src_sqlite(db.file, create = FALSE)
  snpinfo = tbl(my_db, sql(paste0("SELECT * FROM snps WHERE chr='", 
            chr, "' AND pos_Mbp>=", start, " AND pos_Mbp<=", end))) %>%
            collect(n = Inf)
  colnames(snpinfo)[c(1,3)] = c("snp", "pos")
  snpinfo = index_snps(map = map, snpinfo)
  snppr = genoprob_to_snpprob(probs, snpinfo)
  gwas = scan1(genoprobs = snppr, pheno = pheno,
             kinship = K, addcovar = addcovar, intcovar = intcovar, cores = ncores)

  return(list(gwas, snpinfo))

} # assoc_map()


sdp_plot = function(scan1output, snpinfo, thr) {

  map = qtl2plot:::snpinfo_to_map(snpinfo)
  tmp = qtl2plot:::expand_snp_results(scan1output, map, snpinfo)
  snpinfo = cbind(snpinfo, tmp$lod)
  rm(tmp)
  snpinfo.ss = snpinfo[snpinfo[,ncol(snpinfo)] >= thr,]
  sdps = matrix(0, nrow = nrow(snpinfo.ss), ncol = 8, dimnames = 
         list(snpinfo.ss$snp, LETTERS[1:8]))
  unique.sdps = unique(snpinfo.ss$sdp)
  for(j in 1:length(unique.sdps)) {
    wh = which(snpinfo.ss$sdp == unique.sdps[j])
    div = floor(log(unique.sdps[j], 2))
    sdps[wh,div+1] = 1
    rem = unique.sdps[j] - 2^div
    while(rem > 0) {
      div = floor(log(rem, 2))
      sdps[wh,div+1] = 1
      rem = rem - 2^div
    }
  }
  full.sdps = matrix(0, nrow(snpinfo), 8, dimnames = list(snpinfo$snp, LETTERS[1:8]))
  full.sdps[rownames(sdps),] = sdps
  full.sdps = full.sdps[,8:1]
  image(snpinfo$pos, 1:ncol(full.sdps), full.sdps, breaks = c(-0.5, 0.5, 1.5), col = c("white", "black"), axes = F, ann = F)
  wh = which(colSums(full.sdps) > 0)
  for(j in wh) {
    wh2 = which(full.sdps[,j] > 0)
    x = matrix(rep(snpinfo$pos[wh2], each = 2), ncol = 2, byrow = T)
    x = cbind(x, NA)
    y = matrix(rep(c(j-0.5, j+0.5), length(wh2)), ncol = 2, byrow = T)
    y = cbind(y, NA)
    lines(as.vector(t(x)), as.vector(t(y)), lwd = 2)
  }
  abline(h = 1:8 + 0.5, col = "grey80")
  box()
  mtext(side = 2, line = 0.2, at = 8:1, text = names(CCcolors), las = 2)

} # sdp_plot()

get_genes = function(chr, start, end) {

  genes = ensembl[ensembl$type == "gene" & seqnames(ensembl) == chr & 
                  end(ensembl) >= start  * 1e6 & start(ensembl) < end * 1e6]
  return(data.frame(chr = runValue(seqnames(genes)), start = start(genes), 
         stop = end(genes), strand = strand(genes), Name = genes$gene_name))

} # get_genes()

assoc_plot = function(assoc, snpinfo, map, chr, start, end) {

  genes = get_genes(chr, start, end)

  layout(matrix(1:3, 3, 1), heights = c(0.2, 0.3, 0.5))
  par(plt = c(0.08, 0.99, 0, 0.9))
  sdp_plot(scan1output = assoc, snpinfo = snpinfo, thr = max(assoc[,1]) - 1)
  par(plt = c(0.08, 0.99, 0, 1.0))
  plot_snpasso(scan1output = assoc, snpinfo = snpinfo, drop.hilit = 1, 
               xaxt = "n")
  par(plt = c(0.08, 0.99, 0.12, 1))
  plot_genes(genes = genes, xlim = c(start, end), colors = "black")

} # assoc_plot()

gene_snp_intersect = function(genes, scan1output, snpinfo, thr) {

  local.map  = qtl2plot:::snpinfo_to_map(snpinfo)
  local.snps = qtl2plot:::expand_snp_results(scan1output, local.map, snpinfo)
  keep = which(local.snps$lod >= thr)
  local.snps = list(lod = local.snps$lod[keep], map = local.snps$map[keep])
  snpinfo = snpinfo[keep,]
  snp.gr = GRanges(seqnames = snpinfo$chr, ranges = IRanges(start = snpinfo$pos * 1e6, width = 1))
  mcols(snp.gr) = snpinfo
  exons = genes[genes$type %in% "exon"]
  genes.int = subsetByOverlaps(exons, snp.gr)
  snps.int  = subsetByOverlaps(snp.gr, exons)

  snps.int = as.data.frame(snps.int)
  snps.int$symbol = genes.int$gene_name[match(snps.int$ensembl_gene, genes.int$gene_id)]
  snps.int = snps.int[,c(6:9, 11, 16, 12)]
  snps.int = snps.int[snps.int$csq != "intergenic_variant",]
  
  return(snps.int)

} # gene_snp_intersect()
