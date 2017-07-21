#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=3:00:00
module load R/3.3.2

# For each chromosome, harvest the maximum peak per analyte.

cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts

# Arguments.
# INPUTFILE: full path to aggregated QTL LOD file (as *.rds).
# OUTPUTPREFIX: full path to output file prefix. We append .csv and .rds to it.

# JAX normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_liver_metabolites_jax.Rout


# JAX normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_liver_metabolites_jax_sex_gen.Rout


# U. Wisc. normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_liver_metabolites_uwisc.Rout


# U. Wisc. normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_liver_metabolites_uwisc_sex_gen.Rout


# JAX normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_liver_lipids_jax.Rout


# U. Wisc. normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_liver_lipids_uwisc.Rout


# JAX normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_plasma_lipids_jax.Rout


# U. Wisc. normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTFILE=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm_all_qtl.rds
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm
R --no-save --args ${INPUTFILE} ${OUTPUTPREFIX} < harvest_max_qtl.R > harvest_QTL_plasma_lipids_uwisc.Rout

