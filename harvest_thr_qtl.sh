#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=3:59:00
module load R/3.4.1

# For each chromosome, harvest the maximum peak above a given threshold for each analyte.
# Run this AFTER the harvest_max_qtl.R script because it creates the aggregated QTL file.

cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

THR=6

# JAX normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_liver_metabolites_jax.Rout


# JAX normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_liver_metabolites_jax_sex_gen.Rout


# U. Wisc. normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_liver_metabolites_uwisc.Rout


# U. Wisc. normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_liver_metabolites_uwisc_sex_gen.Rout


# JAX normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_liver_lipids_jax.Rout


# U. Wisc. normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_liver_lipids_uwisc.Rout


# JAX normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_plasma_lipids_jax.Rout


# U. Wisc. normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_plasma_lipids_uwisc.Rout


# JAX normalized plasma metabolites. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/metabolites/plasma_metabolites_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/metabolites/plasma_metabolites
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_plasma_metabolites.Rout


# JAX normalized cecum lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/lipids_norm_jax/cecum_lipids_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/lipids_norm_jax/cecum_lipids_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_cecum_lipids_jax.Rout


# JAX normalized cecum  metabolites. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/metabolites/cecum_metabolites_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/metabolites/cecum_metabolites
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_cecum_metabolites.Rout


# JAX normalized islet proteins. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Islet/proteins/islet_proteins_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Islet/proteins/islet_proteins_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} ${THR} < harvest_thr_qtl.R > harvest_qtl_islet_proteins_jax.Rout
