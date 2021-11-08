ref_fasta=$1
chrom_sizes=$2
flank_size=$3
stride=$4
out_prefix="genomewide_gc_hg38_stride_"$stride"_flank_size_"$flank_size

python get_genomewide_gc_bins.py --ref_fasta $ref_fasta \
		  --chrom_sizes $chrom_sizes \
		  --out_prefix $out_prefix \
		  --region_size $flank_size \
		  --stride $stride \
		  --output_format tsv


