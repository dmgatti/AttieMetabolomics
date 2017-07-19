#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=1:00:00
module load R/3.3.2
cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts

R --no-save --args $INDIR $OUTPREFIX $CHUNKSIZE $CHUNKNUM $MAXCOL < qtl2_scan_engine.R > qtl2_scan_engine_${CHUNKNUM}.Rout