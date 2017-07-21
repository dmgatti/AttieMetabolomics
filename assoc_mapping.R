########################
# Association mapping. #
########################
# NOTE: This still a bit messy. Improvements are coming...
# NOTE2: You need /data/gac/resource/CCsnps/ccfoundersnps.sqlite, which is 
#        a database of Sanger SNPs for the 8 DO founders.

# Arguments:
# probs: genoprobs object in qtl2 format.
# pheno: matrix of phenotypes, samples in rows, phenotypes in columns.
# idx: column index of the phenotype that you want to map.
# addcovar: covariates matrix as used in scan1().
# intcovar: covariate to interact with QTL.
# k: list of kinship matrices.
# markers: data.frame containing 4 columns: marker, chr, bp, cM.
# chr: Chromosome to map on.
# start: start position for mapping in Mb.
# end: end position for mapping in Mb.
# ncores: number of cores to use in mappin.
# db.file: Location of the mySQL database containing the Sanger SNPs.
assoc_mapping = function(probs, pheno, idx, addcovar, intcovar = NULL, k, markers, 
                chr, start, end, ncores, 
                db.file = "/data/gac/resource/CCsnps/ccfoundersnps.sqlite") {

  # Subset probs and K to keep only the current chromosome.
  probs = probs[,chr]
  k     = k[[chr]]

  # Split up markers into a vector of map positions.
  map = split(markers[,3] * 1e-6, markers[,2])
  nm  = split(markers[,1], markers[,2])
  map = mapply(function(x, y) { names(x) = y;x }, map, nm)
  map = map[order(as.numeric(names(map)))]

  # Extract SNPs from the database
  my_db = src_sqlite(db.file, create = FALSE)
  snpinfo = tbl(my_db, sql(paste0("SELECT * FROM snps WHERE chr='", 
            chr, "' AND pos_Mbp>=", start, " AND pos_Mbp<=", end))) %>%
            collect(n = Inf)

  # Names have to be replaced for future methods
  colnames(snpinfo)[c(1,3)] = c("snp", "pos")

  # Index groups of similar SNPs.
  snpinfo = index_snps(map = map, snpinfo)

  # Find which phenotype data actually exist
  sel = !is.na(pheno[,idx])

  # Convert genoprobs to snpprobs.
  snppr = genoprob_to_snpprob(probs[sel,], snpinfo)
  
  # Scan1.
  assoc = scan1(pheno = pheno[sel,idx, drop = FALSE], kinship = k[sel,sel],
          genoprobs = snppr, addcovar = addcovar[sel,], 
          intcovar = addcovar[sel,intcovar], cores = ncores)

  # Return the scan data.
  return(list(assoc, snpinfo))

} # assoc_mapping()
