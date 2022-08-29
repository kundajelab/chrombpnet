input_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/sorted_merged.bam
output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/



#samtools view -H $input_bam > $output_dir"/replicates/pr_header.sam"

##Split merged treatment BAM
#nlines=$(samtools view $input_bam | wc -l ) # Number of reads in the BAM file
#nlines=$(( (nlines + 1) / 2 )) # half that number
#samtools view $input_bam | shuf - | split -d -l ${nlines} - "$output_dir/replicates/pr_input" # This will shuffle the lines in the file and split in two 
#cat $output_dir/replicates/pr_header.sam $output_dir/replicates/pr_input00 | samtools view -bS - > $output_dir/replicates/pr_input00.bam
#cat $output_dir/replicates/pr_header.sam $output_dir/replicates/pr_input01 | samtools view -bS - > $output_dir/replicates/pr_input01.bam



blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes


in_bam=$output_dir/replicates/pr_input00.bam
samtools sort  $output_dir/replicates/pr_input00.bam -o $output_dir/replicates/pr_input00.sorted.bam
in_bam=$output_dir/replicates/pr_input00.sorted.bam

bash step2_make_bigwigs_from_bams.sh $in_bam $output_dir/replicates/pr1 ATAC_PE $ref_fasta $chrom_sizes 

#in_bam=$output_dir/replicates/pr_input01.bam

in_bam=$output_dir/replicates/pr_input01.bam
samtools sort  $output_dir/replicates/pr_input01.bam -o $output_dir/replicates/pr_input01.sorted.bam
in_bam=$output_dir/replicates/pr_input01.sorted.bam

bash step2_make_bigwigs_from_bams.sh $in_bam $output_dir/replicates/pr2 ATAC_PE $ref_fasta $chrom_sizes 
