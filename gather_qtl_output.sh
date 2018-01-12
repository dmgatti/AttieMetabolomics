#!/bin/bash -l
#PBS -q short -l nodes=cadillac022:ppn=10,walltime=3:59:00
module load R/3.4.1

# For each chromosome, harvest the maximum peak per analyte.

cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

# JAX normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/liver_metabolites_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_metabolites_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_jax.Rout


# JAX normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/liver_metabolites_jax_norm_sex_gen
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_metabolites_norm_jax_sex_gen/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_jax_sex_gen.Rout


# U. Wisc. normalized liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc/liver_metabolites_uwisc_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_metabolites_norm_uwisc/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_uwisc.Rout


# U. Wisc. normalized liver metabolites with sex & gen as covariates. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_uwisc_sex_gen/liver_metabolites_uwisc_norm_sex_gen
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_metabolites_norm_uwisc_sex_gen/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_uwisc_sex_gen.Rout


# JAX normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_lipids_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_lipids_jax.Rout


# U. Wisc. normalized liver lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_uwisc/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_uwisc/liver_lipids_uwisc_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_lipids_norm_uwisc/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_lipids_uwisc.Rout


# JAX normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_jax/plasma_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/plasma_lipids_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_lipids_jax.Rout


# U. Wisc. normalized plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_uwisc/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/lipids_norm_uwisc/plasma_lipids_uwisc_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/plasma_lipids_norm_uwisc/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_lipids_uwisc.Rout


# JAX normalized plasma metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/metabolites/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma/metabolites/plasma_metabolites
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/plasma_metabolites/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_metabolites_jax.Rout


# JAX normalized cecum lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/lipids_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/lipids_norm_jax/cecum_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/cecum_lipids_norm_jax/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_cecum_lipids_jax.Rout


# JAX normalized cecum metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/metabolites/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum/metabolites/cecum_metabolites
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/cecum_metabolites/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_cecum_metabolites_jax.Rout


# JAX normalized islet proteins. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Islet/proteins/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Islet/proteins/islet_proteins_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/islet_proteins/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_islet_proteins_jax.Rout

###############
# RankZ phenotypes.

# Cecum lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum_rankZ/lipids/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum_rankZ/lipids/cecum_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/Cecum_rankZ/lipids/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_cecum_lipids_jax.Rout

# Cecum metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum_rankZ/metabolites/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Cecum_rankZ/metabolites/cecum_metabolites_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/Cecum_rankZ/metabolites/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_cecum_metabolites_jax.Rout

# Liver lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver_rankZ/lipids_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver_rankZ/lipids_norm_jax/liver_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/Liver_rankZ/lipids/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_lipids_jax.Rout

# Liver metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver_rankZ/metabolites_norm_jax/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver_rankZ/metabolites_norm_jax/liver_metabolites_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/Liver_rankZ/metabolites/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_liver_metabolites_jax.Rout

# Plasma lipids. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma_rankZ/lipids/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma_rankZ/lipids/plasma_lipids_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/Plasma_rankZ/lipids/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_lipids_jax.Rout

# Plasma metabolites. (NOTE: Place a / at the end of the paths)
INPUTDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma_rankZ/metabolites/
OUTPUTPREFIX=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma_rankZ/metabolites/plasma_metabolites_jax_norm
FIGDIR=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/Plasma_rankZ/metabolites/
R --no-save --args ${INPUTDIR} ${OUTPUTPREFIX} ${FIGDIR} < gather_qtl_output.R > gather_qtl_plasma_metabolites_jax.Rout




