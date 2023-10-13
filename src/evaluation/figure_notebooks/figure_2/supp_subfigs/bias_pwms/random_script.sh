
genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes

python random_build_pwm_from_bigwig.py  -i /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw -g $genome -o random.pwm -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 

