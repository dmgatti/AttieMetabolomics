cd /hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/scripts

INPUT=/hpcdata/gac/derived/Attie_DO_Liver_Metabolomics/qtl2_input/attie_liver_metabolites_qtl2_input_sex_gen.Rdata
OUTPUT=/hpcdata/gac/projects/Attie_DO_Liver_Metabolomics/QTL/Liver/metabolites_norm_jax_sex_gen/

# Set the maximum number in the loop to ceiling(MAXCOL/10).
for i in {16..31}
do
  qsub -v INDIR=${INPUT},OUTPREFIX=${OUTPUT},CHUNKSIZE=10,CHUNKNUM=${i},MAXCOL=308 qtl2_scan_engine.sh
done
