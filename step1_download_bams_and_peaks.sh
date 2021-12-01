# make data directory
mkdir data

# download bam
wget https://www.encodeproject.org/files/ENCFF415FEC/@@download/ENCFF415FEC.bam -O data/rep1.bam
wget https://www.encodeproject.org/files/ENCFF646NWY/@@download/ENCFF646NWY.bam -O data/rep2.bam
samtools merge -f data/merged.bam data/rep1.bam data/rep2.bam
samtools index data/merged.bam

# download overlap peaks
wget https://www.encodeproject.org/files/ENCFF478ZAW/@@download/ENCFF478ZAW.bed.gz -O data/overlap.bed.gz

# download idr peaks
wget https://www.encodeproject.org/files/ENCFF945SYZ/@@download/ENCFF945SYZ.bed.gz  -O  data/idr.bed.gz

#download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O data/hg38.fa.gz
gunzip data/hg38.fa.gz

# download reference chromsome sizes
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O data/hg38.chrom.sizes

# download blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O data/blacklist.bed.gz