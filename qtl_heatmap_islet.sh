#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=12:00:00
cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

module load R/3.4.1

# Arguments:
# INFILE:  full path to the LOD matrix for all phenotypes, stored as a *.rds file.
# OUTFILE: full path to the figure file, saved as PNG.
# LODTHR:  LOD threshold above which LOD scores will be truncated.

LODTHR=8


##########
# Islet proteins, JAX normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Islet/proteins/islet_proteins_jax_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/islet_proteins/islet_proteins_jax_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_islet_proteins_jax_norm.Rout
