
file=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/interpret_uncorrected_model_05.10.2022/GM12878
chrom_sizes=/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes
python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1


file=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/interpret_DNASE_SE_03.06.2022_simplebias/GM12878
chrom_sizes=/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes
python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1


#file=/mnt/lab_data2/anusri/chrombpnet/results/hint_atac/DNASE_SE/GM12878/DNASE_SE_11.28.2022_hint_atac/hint_atac_model/interpret/GM12878
#chrom_sizes=/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes
#python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1


#file=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/interpret/GM12878
#chrom_sizes=/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes
#python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

#file1=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/DNASE_SE/GM12878/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/BIAS/GM12878
#chrom_sizes=/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes
#python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file1.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file1.profile.bw -s $file1.profile.stat -t 1





