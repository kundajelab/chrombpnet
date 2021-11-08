#get the inverse intersection of idr peak file and all gc genome bins
task=$1
idr=$2
genomewide_gc=$3
chrom_sizes=$4
flank_size=$5
blacklist=$6

echo $blacklist
bedtools slop -i $blacklist -g $chrom_sizes -b $flank_size > $task/blacklist_slop1024.bed
bedtools slop -i $idr -g $chrom_sizes -b $flank_size | bedtools intersect -v -a $genomewide_gc -b stdin $task/blacklist_slop1024.bed > $task/candidate.negatives.tsv
rm $task/blacklist_slop1024.bed
