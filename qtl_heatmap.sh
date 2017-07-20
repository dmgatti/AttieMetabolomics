#!/bin/bash -l
cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts

# Arguments:
# INFILE:  full path to the LOD matrix for all phenotypes, stored as a *.rds file.
# OUTFILE: full path to the figure file, saved as PNG.
# LODTHR:  LOD threshold above which LOD scores will be truncated.

##########
# Liver lipids, JAX normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_lipids_norm_jax/liver_lipids_jax_qtl_heatmap.png
LODTHR=8

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_liver_lipids_jax_norm.Rout

##########
# Liver lipids, U. Wisc. normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_lipids_norm_uwisc/liver_lipids_uwisc_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_liver_lipids_uwisc_norm.Rout

##########
# Liver metabolites, JAX normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_jax/liver_metabolites_jax_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_liver_metabolites_jax_norm.Rout

##########
# Liver metabolites, JAX normalized, sex & gen.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_jax_sex_gen/liver_metabolites_jax_sex_gen_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_liver_metabolites_sex_gen_jax_norm.Rout

##########
# Liver metabolites,  U. Wisc. normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_uwisc/liver_metabolites_uwisc_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_liver_metabolites_uwisc_norm.Rout

##########
# Liver metabolites,  U. Wisc. normalized, sex & gen.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_sex_gen_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_liver_metabolites_uwisc_sex_gen_norm.Rout

##########
# Plasma lipids, JAX normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/plasma_lipids_norm_jax/plasma_lipids_jax_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_plasma_lipids_jax_norm.Rout

##########
# Plasma lipids, U. Wisc. normalized, sex, gen & batch.
INFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm_all_qtl.rds
OUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/palsma_lipids_norm_uwisc/plasma_lipids_uwisc_qtl_heatmap.png

R --no-save --args ${INFILE} ${OUTFILE} ${LODTHR} < qtl_heatmap.R > qtl_heatmap_plasma_lipids_uwisc_norm.Rout

