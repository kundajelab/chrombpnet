#get the inverse intersection of idr peak file and all gc genome bins
task=$1
idr=$2
genomewide_gc=$3
chrom_sizes=$4
flank_size=$5

bedtools slop -i $idr -g $chrom_sizes -b $flank_size | bedtools intersect -v -a $genomewide_gc -b stdin > $task/candidate.negatives.tsv
