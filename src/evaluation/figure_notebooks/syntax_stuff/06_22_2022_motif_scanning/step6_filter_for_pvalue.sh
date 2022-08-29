#dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/moods/profile/
dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/moods/counts/


awk -F'\t' '{if($13<0.001 || $13>0.999) {print $0}}' $dir/moods_modisco_hits_with_pvalues.bed | bedtools sort -i stdin | bgzip >  $dir/moods_modisco_hits_with_pvalues_0.001_sig.bed.gz
awk -F'\t' '{if($13<0.01 || $13>0.99) {print $0}}' $dir/moods_modisco_hits_with_pvalues.bed | bedtools sort -i stdin | bgzip >  $dir/moods_modisco_hits_with_pvalues_0.01_sig.bed.gz



