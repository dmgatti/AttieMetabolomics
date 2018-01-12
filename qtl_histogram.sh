#!/bin/bash -l
module load R/3.4.1

# Arguments:
# input.file: full path to the *.rds qtl summary file.
# output.file: full path to the output figure file as a PNG.
# thr: LOD threshold to use when selecting QTL peaks.

cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts
BASEDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/


#### RankZ phenotypes ####


##########
# Liver lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Liver_rankZ/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/Liver_rankZ/lipids/liver_lipids_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_liver_lipids.Rout

##########
# Liver metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Liver_rankZ/metabolites_norm_jax/liver_metabolites_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/Liver_rankZ/metabolites/liver_metabolites_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_liver_metabolites.Rout

##########
# Plasma lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Plasma_rankZ/lipids/plasma_lipids_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/Plasma_rankZ/lipids/plasma_lipids_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_plasma_lipids.Rout

##########
# Plasma metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Plasma/metabolites/plasma_metabolites_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/Plasma_rankZ/metabolites/plasma_metabolites_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_plasma_metabolites.Rout

##########
# Cecum lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Cecum_rankZ/lipids/cecum_lipids_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/Cecum_rankZ/lipids/cecum_lipids_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_cecum_lipids.Rout

##########
# Cecum metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Cecum_rankZ/metabolites/cecum_metabolites_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/Cecum_rankZ/metabolites/cecum_metabolites_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_cecum_metabolites.Rout





#### Non-rankZ pheotypes ####

##########
# Liver lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/liver_lipids_norm_jax/liver_lipids_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_liver_lipids.Rout

##########
# Liver metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/liver_metabolites_norm_jax/liver_metabolites_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_liver_metabolites.Rout

##########
# Plasma lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/plasma_lipids_norm_jax/plasma_lipids_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_plasma_lipids.Rout

##########
# Plasma metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Plasma/metabolites/plasma_metabolites_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/plasma_metabolites/plasma_metabolites_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_plasma_metabolites.Rout

##########
# Cecum lipids: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Cecum/lipids_norm_jax/cecum_lipids_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/cecum_lipids_norm_jax/cecum_lipids_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_cecum_lipids.Rout

##########
# Cecum metabolites: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Cecum/metabolites/cecum_metabolites_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/cecum_metabolites/cecum_metabolites_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_plasma_metabolites.Rout


##########
# Islet Proteins: JAX: sex, gen & batch
INFILE=${BASEDIR}QTL/Islet/proteins/islet_proteins_jax_norm_qtl_summary_thresh_6.csv
OUTFILE=${BASEDIR}figures/QTL/islet_proteins/islet_proteins_jax_norm_qtl_histogram.png

R --no-save --args ${INFILE} ${OUTFILE} < qtl_histogram.R > qtl_histogram_islet_proteins.Rout

