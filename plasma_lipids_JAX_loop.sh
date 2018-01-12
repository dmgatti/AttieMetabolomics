cd /hpcdata/gac/projects/Attie_DO_Metabolomics/scripts

INPUT=/hpcdata/gac/derived/Attie_DO_Metabolomics/qtl2_input/attie_plasma_lipids_qtl2_input.Rdata
OUTPUT=/hpcdata/gac/projects/Attie_DO_Metabolomics/QTL/Plasma_rankZ/lipids/
SIZE=20
RZ=TRUE

# Set the maximum number in the loop to ceiling(MAXCOL/SIZE).
for i in {1..88}
do
  echo $i
  qsub -v INDIR=${INPUT},OUTPREFIX=${OUTPUT},CHUNKSIZE=${SIZE},CHUNKNUM=${i},RANKZ=${RZ} qtl2_scan_engine.sh
  sleep 30s
done
