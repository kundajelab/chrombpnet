ref_fasta=$1
chrom_sizes=$2
inputlen=$3
stride=$4
out_prefix="genomewide_gc_hg38_stride_"$stride"_inputlen_"$inputlen

python get_genomewide_gc_bins.py --ref_fasta $ref_fasta \
		  --chrom_sizes $chrom_sizes \
		  --out_prefix $out_prefix \
		  --region_size $inputlen \
		  --stride $stride \
		  --output_format tsv


