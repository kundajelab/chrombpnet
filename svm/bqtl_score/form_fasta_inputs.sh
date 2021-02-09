#note: the Frazer bQTL files are hg19, 1-indexed 
#run in gm12878, fold 0
for tf in junD nfkB oct11 oct1 pu1 stat1
do
    #junD.formatted.csv  nfkB.formatted.csv	oct11.formatted.csv  oct1.formatted.csv  pu1.formatted.csv  stat1.formatted.csv
    echo $tf
    python form_fasta_inputs.py --snp_file src_formatted/$tf.formatted.csv \
	   --chrom_col Chr \
	   --pos_col position \
	   --ref_col POSTallele \
	   --alt_col ALTallele \
	   --ref_fasta /mnt/data/male.hg19.fa \
	   --offset_1 \
	   --rsid_col rsid \
	   --out_prefix gm12878_gkm_preds_bqtl_$tf &
done
