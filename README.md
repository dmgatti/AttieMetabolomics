Metabolite and lipid data from Diversity Outbred mice.

Liver metabolites: 308 phenotypes.
Liver lipids: 1354 phenotypes.
Plasma lipids: 1767 phenotypes.

***

## QTL Mapping Pipeline

1. Format and normalize the raw analyte data.
  + Modify file 'attie_liver_metabolites_normalize.R'.
  + If there is no missing data, you don't need to impute data, just batch normalize using ComBat.
  + Name the first column containing the mouse IDs "Mouse.ID" (case-sensitive).
  + Save a \*.rds file containing the normalized analytes.  
2. Gather the QTL input data into a single file.
  + Modify file 'gather_qtl_input_data.R'.
  + Read in the analyte file produced in setp 1 and the genoprobs, located on the [JAX FTP site](ftp://ftp.jax.org/dgatti/_forKarl/attie_do_genoprobs_20170522.rds).
  + This script will create a single compressed R binary file (\*.Rdata) containing objects called:
    + pheno: data.frame containing the normalized phenotypes and covariates.
    + pheno.rz: data.frame containing the Z-score transformed phenotypes and covariates.
    + pheno.descr: a small phenotype dictionary.
    + genoprobs: Haplotype probabilities for all mice in qtl2 format.
    + K: list of kinship matrices.
    + map: data.frame containing marker information.
 3. Map the analytes.
   + Modify file 'qtl2_scan_engine.R'.
   + The script is set up to run 1000's of analytes in chunks on the cluster.
   + You can run all analytes together as long as they all use the same covariates.
   + The output is a \*.rds file containing LOD scores for all analytes (num_markers X num_analytes).
 4. Collect the LOD scores for all chunks into one file.
   + Modify file 'gather_qtl_output.R'.
   + This file will gather the chunked QTL files an write them out to a single \*.rds file.
   + It will also create a PDF containing all QTL plots.
 5. Harvest the QTL peaks above your desired threshold.
   + Modify 'harvest_thr_qtl.R'.
   + The input file is the \*.rds file containing all LOD scores created in step 4.
   + Set a LOD threshold.
   + The output will be a \*.csv file containing the highest peak on each chromosome above the LOD threshold.
 6. Create QTL heatmap and histogram.
   + Modify 'qtl_heatmap.R' and 'qtl_histogram.R'.
   + The input for the QTL heatmap is the \*.rds file created in step 4 containing all LOD scores at all markers.
   + The output is a large PNG.
   + The input for the QTL histogram is the \*.csv file created in step 5 containing the harvested QTL peaks.
   + The output is a PNG that plots the chromosomes that have peaks within 2 Mb.
 7. Calculate founder allele effects and association mapping LODs at selected peaks.
  + Modify 'qtl2_coef_assoc_engine.R'.
  + The input is the \*.rds file created in step 4 and the harvested QTL file created in step 5.
  + The output will be one founder allele effects plots and one association mapping plot for each peak in the harvested QTL file.
  
  
