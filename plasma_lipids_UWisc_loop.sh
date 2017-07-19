cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts

INPUT=/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/qtl2_input/attie_plasma_lipids_qtl2_UWisc_input.Rdata
OUTPUT=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Plasma/lipids_norm_uwisc/

# Set the maximum number in the loop to ceiling(MAXCOL/10).
for i in {1..177}
do
  qsub -v INDIR=${INPUT},OUTPREFIX=${OUTPUT},CHUNKSIZE=10,CHUNKNUM=${i},MAXCOL=1767 qtl2_scan_engine.sh
done
