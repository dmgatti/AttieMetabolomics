# Get the lods for a set of gene, metabolite or lipid QTL.
# path: character string with full path to LOD file.
# analytes: character vector containing a set of ensembl IDs or analytes.
get_lods = function(path, analytes) {

  # Get the full LOD file.
  lod = readRDS(path)

  # Check for analytes that aren't in the LOD file.
  wh = which(!analytes %in% colnames(lod))
  if(length(wh) > 0) {

    print(paste("Some analytes not found in LOD file:", 
          paste(analytes[wh], collapse = ", ")))

  } # if(length(wh) > 0)

  return(lod[,analytes])

} # get_lods()


# Get the coefficients fora given chromosome and set of genes or analytes.
# path: character string with full path to coef directory.
# analytes: data.frame with 2 columns: 
#           analyte: character vector containing a set of ensembl IDs or analytes,
#           chr: character vector containing chromosome IDs.
get_coefs = function(path, analytes) {

  coef.pattern = paste0(analytes$analyte, "_coef_chr", analytes$chr)

  filename = dir(path = path, pattern = coef.pattern[1], full.names = TRUE)
  x = readRDS(filename)

  mat = matrix(0, nrow = 8 * nrow(x), ncol = length(coef.pattern), dimnames = 
        list(paste(rownames(x), rep(LETTERS[1:8], each = nrow(x)), sep = "_"), 
        analytes$analyte))

  for(i in 1:length(coef.pattern)) {

    filename = dir(path = path, pattern = coef.pattern[i], full.names = TRUE)
    if(length(filename) > 0) {

      x = readRDS(filename)
      mat[,i] = as.vector(x[,1:8])

    } else {

      print(paste(coef.pattern[i], "not found."))

    } # else

  } # for(i)

  return(mat)

} # get_coefs()
