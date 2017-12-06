#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=3:59:00
module load R/3.4.1
cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

R CMD BATCH --no-save gather_qtl_input_data.R

