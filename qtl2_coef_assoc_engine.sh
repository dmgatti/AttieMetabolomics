#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=10,walltime=3:59:00
module load R/3.4.1
cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

# DF: The full path to the qtl2 data input file. *.Rdata
# LF: The full path to the QTL LOD file. *.rds
# QF: The full path to the QTL summary file. *.csv
# OUTDIR: Output directory for data.
# FIGDIR: Output directory for figures.
# CHUNKSIZE: Number of analytes to map in one chunk.
# CHUNKNUM:  Index of chunk to run.

R --no-save --args $DATAFILE $LODFILE $QTLFILE $OUTDIR $FIGDIR $CHUNKSIZE $CHUNKNUM < qtl2_coef_assoc_engine.R > qtl2_coef_assoc_engine_${CHUNKNUM}.Rout
