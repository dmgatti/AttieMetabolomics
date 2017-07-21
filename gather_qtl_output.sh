#!/bin/bash -l
module load R/3.3.2

# For each chromosome, harvest the maximum peak per analyte.

cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts

# JAX normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_jax.Rout


# JAX normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_jax_sex_gen/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_jax_sex_gen.Rout


# U. Wisc. normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_uwisc/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_uwisc.Rout


# U. Wisc. normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_metabolites_norm_uwisc_sex_gen/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_uwisc_sex_gen.Rout


# JAX normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_lipids_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_lipids_jax.Rout


# U. Wisc. normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_uwisc/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/liver_lipids_norm_uwisc/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_lipids_uwisc.Rout


# JAX normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/plasma_lipids_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_lipids_jax.Rout


# U. Wisc. normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_uwisc/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/figures/QTL/plasma_lipids_norm_uwisc/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_lipids_uwisc.Rout

