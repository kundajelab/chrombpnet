input_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/sorted_merged.bam
output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/

#name sort to achieve random genomic pos#
samtools sort -n --threads 50 -O bam -o $output_dir/replicates/namesorted.bam $input_bam
echo "done sorting" 


samtools view -H $output_dir/replicates/namesorted.bam > $output_dir/replicates/pr1.sam
cp $output_dir/pr1.sam $output_dir/replicates/pr2.sam

echo "made headers"
sam1=$output_dir/replicates/pr1.sam
sam2=$output_dir/replicates/pr2.sam 
samtools view --threads 50  $output_dir/replicates/namesorted.bam | awk -v sam1="$sam1" -v sam2="$sam2" '{if(NR%2){print >> sam1} else {print >> sam2}}'

samtools view -bS $sam1 > $output_dir/replicates/pr1.bam
samtools view -bS $sam2 > $output_dir/replicates/pr2.bam

blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes

in_bam=$output_dir/replicates/pr1.bam
bash step2_make_bigwigs_from_bams.sh $in_bam $output_dir/replicates/pr1 ATAC_PE $ref_fasta $chrom_sizes 

in_bam=$output_dir/replicates/pr2.bam
bash step2_make_bigwigs_from_bams.sh $in_bam $output_dir/replicates/pr2 ATAC_PE $ref_fasta $chrom_sizes 
