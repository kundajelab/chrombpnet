data_dir=$1

# download bam
wget https://www.encodeproject.org/files/ENCFF415FEC/@@download/ENCFF415FEC.bam -O $data_dir/rep1.bam
wget https://www.encodeproject.org/files/ENCFF646NWY/@@download/ENCFF646NWY.bam -O $data_dir/rep2.bam
samtools merge -f $data_dir/merged.bam $data_dir/rep1.bam $data_dir/rep2.bam
samtools index $data_dir/merged.bam

# download overlap peaks
wget https://www.encodeproject.org/files/ENCFF470YYO/@@download/ENCFF470YYO.bed.gz -O $data_dir/overlap.bed.gz

#download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O $data_dir/hg38.fa.gz
gunzip $data_dir/hg38.fa.gz

# download reference chromsome sizes
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O $data_dir/hg38.chrom.sizes

# download blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O $data_dir/blacklist.bed.gz
