#get the inverse intersection of idr peak file and all gc genome bins
task=$1
idr=$2
genomewide_gc=$3
bedtools intersect -v -a $genomewide_gc -b $idr > $task/candidate.negatives.tsv
