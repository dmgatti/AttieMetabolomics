#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=2:00:00
module load R/3.4.1
cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

R --no-save --args $INDIR $OUTPREFIX $CHUNKSIZE $CHUNKNUM $RANKZ < qtl2_scan_engine.R > qtl2_scan_engine_${CHUNKNUM}.Rout
