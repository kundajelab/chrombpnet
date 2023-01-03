output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/NaiveCD4/ATAC_PE_09.21.2022_withgm12878bias/SIGNAL_full/modisco_crop_500_lite_200K/annotations/
peak_bed=/oak/stanford/groups/akundaje/anusri/Mineto-data/NaiveCD4/ATAC_PE_09.21.2022_withgm12878bias/chrombpnet_model/interpret_full/NaiveCD4.interpreted_regions_v2.bed
modisco_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/NaiveCD4/ATAC_PE_09.21.2022_withgm12878bias/SIGNAL_full/modisco_crop_500_lite_200K/
shap_dir=/oak/stanford/groups/akundaje/anusri/Mineto-data/NaiveCD4/ATAC_PE_09.21.2022_withgm12878bias/chrombpnet_model/interpret_full/NaiveCD4
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa 


mkdir  $output_dir
echo $output_dir

if [[ -f $output_dir/equal.spaced.merged.bed ]];then
        echo "peaks file exists"
        wc -l $output_dir/equal.spaced.merged.bed
else
        # remove blacklist regions from peaks and take 1000bp around the summit
        bedtools slop -i /mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz -g /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes -b 1057 > $output_dir/temp.bed
        bedtools intersect -v -a $peak_bed -b $output_dir/temp.bed | awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' > $output_dir/merged_peaks_no_blacklist.w1000.bed

        # merge your peak set
        bedtools sort -i $output_dir/merged_peaks_no_blacklist.w1000.bed | bedtools merge -i stdin   > $output_dir/merged.bed

        # create equal intervals from your merged peaks - this is optional - this will affect how you consider co-occurence of hits not the annotation or the output hit bed file itself.
        python create_equal_width_peaks.py --bed $output_dir/merged.bed --outf $output_dir/equal.spaced.merged.bed --input_len 1000

        peaks_bed=$output_dir/equal.spaced.merged.bed
        wc -l $output_dir/equal.spaced.merged.bed
fi

# update to the new peak file
peak_bed=$output_dir/equal.spaced.merged.bed


bash run_moods_hits_test.sh $modisco_dir $shap_dir $output_dir $peak_bed "mean_norm"

