in_bam=$1
interm=$2
chrom_sizes=$3

samtools view -b -@50 -F780 -q30  $in_bam | bedtools bamtobed -i stdin | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2} else if ($6 == "-") {$3 = $3 + 1} print $0}' | bedtools genomecov -5 -i stdin -g $chrom_sizes | bedGraphToBigWig stdin $chrom_sizes $interm/shifted.sorted.bam.bpnet.unstranded.bw

