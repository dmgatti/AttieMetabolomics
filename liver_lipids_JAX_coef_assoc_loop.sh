cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

# The full path to the qtl2 data input file.
DF=/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_liver_lipids_qtl2_input.Rdata
# The full path to the QTL LOD file.
LF=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_all_qtl.rds
# The full path to the QTL summary file.
QF=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/liver_lipids_jax_norm_qtl_summary_thresh_6.csv
# Output directory for data.
OD=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/lipids_norm_jax/
# Output directory for figures.
FD=/hpcdata/gac/projects/Attie_DO_Metabolomics/figures/QTL/liver_lipids_norm_jax/
# Number of analytes to map in one chunk.
SIZE=20


# Set the maximum number in the loop to ceiling(number of rows in QTL summary/SIZE).
# 1774 rows in liver_lipids_jax_norm_qtl_summary_thresh_6.csv
for i in {1..89}
do
  qsub -v DATAFILE=${DF},LODFILE=${LF},QTLFILE=${QF},OUTDIR=${OD},FIGDIR=${FD},CHUNKSIZE=${SIZE},CHUNKNUM=${i} qtl2_coef_assoc_engine.sh
  sleep 30s
done
