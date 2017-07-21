#!/bin/bash -l
module load R/3.3.2

# Arguments:
# input.file: full path to the *.rds qtl summary file.
# output.file: full path to the output figure file as a PNG.
# thr: LOD threshold to use when selecting QTL peaks.

cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts
BASEDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/

THR=6

##########
# Liver lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_all_qtl.rds
OUTFILE=${BASEDIR}figures/QTL/liver_lipids_norm_jax/liver_lipids_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} ${THR} < qtl_histogram.R > qtl_histogram_liver_lipids.Rout

##########
# Liver metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_all_qtl.rds
OUTFILE=${BASEDIR}figures/QTL/liver_metabolites_norm_jax/liver_metabolites_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} ${THR} < qtl_histogram.R > qtl_histogram_liver_metabolites.Rout

##########
# Plasma metabolites: JAX: sex, gen & batc
INFILE=${BASEDIR}QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm_all_qtl.rds
OUTFILE=${BASEDIR}figures/QTL/plasma_lipids_norm_jax/plasma_lipids_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} ${THR} < qtl_histogram.R > qtl_histogram_plasma_lipids.Rout

