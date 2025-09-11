#/oak/stanford/groups/akundaje/shared/nakedDNA/rep1_S1_001.st.rmdup.flt.bam

#out=/oak/stanford/groups/akundaje/shared/nakedDNA/
#samtools merge -f merged_unsorted.bam $out/rep1_S1_001.st.rmdup.flt.bam  $out/rep2_S2_001.st.rmdup.flt.bam $out/rep3_S3_001.st.rmdup.flt.bam $out/rep4_S4_001.st.rmdup.flt.bam $out/rep5_S5_001.st.rmdup.flt.bam
#samtools sort -@4 merged_unsorted.bam -o merged.bam
#samtools index merged.bam

#samtools view -c  merged.bam

data_dir="/mnt/lab_data2/anusri/print_analysis/data"
data_type="ATAC_PE"
ref_fasta=reference/hg38.genome.fa
chrom_sizes=reference/chrom.sizes
in_bam=/mnt/lab_data2/anusri/print_analysis/merged.bam

bash step2_make_bigwigs_from_bams.sh $in_bam $data_dir"/" $data_type $ref_fasta $chrom_sizes 
