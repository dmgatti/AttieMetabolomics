cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

INPUT=/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_liver_metabolites_qtl2_input.Rdata
OUTPUT=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Liver_rankZ/metabolites_norm_jax/
SIZE=20
RZ=TRUE

# Set the maximum number in the loop to ceiling(MAXCOL/SIZE).
for i in {1..15}
do
  echo $i
  qsub -v INDIR=${INPUT},OUTPREFIX=${OUTPUT},CHUNKSIZE=${SIZE},CHUNKNUM=${i},RANKZ=${RZ} qtl2_scan_engine.sh
  sleep 15s
done
