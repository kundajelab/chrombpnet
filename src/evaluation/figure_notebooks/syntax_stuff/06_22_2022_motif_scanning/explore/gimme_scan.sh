peak_bed=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/nautilus_runs_jun16/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/interpret/merged.HEPG2.interpreted_regions.bed
#peak_bed=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/data/peaks_no_blacklist.bed
out_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/
#echo $out_dir

#sort  --parallel=10 -S 30%  -k 8gr,8gr $peak_bed > $out_dir/peaks_sorted.bed
#head -n 500 $out_dir/peaks_sorted.bed | awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' > $out_dir/top500.w1000.bed

#awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' $peak_bed > $out_dir/dnase_peaks_no_blacklist.w1000.bed


## archtype motifs

#bedtools intersect -a /mnt/data/vierstra_motifs/hg38.archetype_motifs.v1.0.bed.gz -b $peak_bed -wa | uniq > $out_dir/merged.archetype.intersect.hepg2.bed
#wc -l $out_dir/merged.archetype.intersect.hepg2.bed


#wc -l $out_dir/dnase_peaks_no_blacklist.w1000.bed
#bed_file=$out_dir/top500.w1000.bed

bed_file=$out_dir/dnase_peaks_no_blacklist.w1000.bed
fasta=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
#out_file=$out_dir/gimme_modisco_hits_hococomo.bed

out_file=$out_dir/gimme_modisco_hits_meme.bed
#ppm_file=$out_dir/counts.ppms.txt
ppm_file=/oak/stanford/groups/akundaje/soumyak/motifs/motifs.meme.txt

#gimme scan $bed_file -p $ppm_file -g $fasta -n 100 -N 10 -f 0.20 -z -b > $out_file &
#gimme scan $bed_file -p $ppm_file -g $fasta -n 100 -N 10 -f 0.20 -z -b > $out_file &

#gimme scan $bed_file -p HOCOCOMOv11_HUMAN -g $fasta -n 100 -N 10 -f 0.50 -z -b > $out_file &





