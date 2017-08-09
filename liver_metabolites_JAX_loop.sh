cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

INPUT=/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_liver_metabolites_qtl2_input.Rdata
OUTPUT=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver/metabolites_norm_jax/
SIZE=20

# Set the maximum number in the loop to ceiling(MAXCOL/SIZE).
for i in {1..16}
do
  qsub -v INDIR=${INPUT},OUTPREFIX=${OUTPUT},CHUNKSIZE=${SIZE},CHUNKNUM=${i} qtl2_scan_engine.sh
  sleep 30s
done
