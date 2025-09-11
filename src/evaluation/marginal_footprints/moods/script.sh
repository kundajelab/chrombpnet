#outdir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/GM12878/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/SIGNAL/modisco_crop_500
#modisco convert -i $outdir/modisco_results_allChroms_counts.hdf5 -o $outdir/modiscolite_counts.h5
#modisco meme -i $outdir/modiscolite_counts.h5 -o "PFM" -t "PFM"
#echo "modisco meme -i $outdir/modiscolite_counts.h5 -o "PFM" -t "PFM""

#samtools faidx /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa chr1 > chr1.fa

#outdir=/mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/moods/
#fimo  --thresh 0.001 --no-qvalue --verbosity 5 --max-stored-scores 2000000000 --parse-genomic-coord  -o $outdir/fimo_output PFM chr1.fa

file1=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed
file2=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR000EMT/preprocessing/downloads/peaks.bed.gz 

zcat $file2 | cat - $file1 > concatenated.bed

# Sort the concatenated BED file
bedtools sort -i concatenated.bed > sorted.bed

# Merge the sorted BED file
bedtools merge -i sorted.bed > merged_output.bed






