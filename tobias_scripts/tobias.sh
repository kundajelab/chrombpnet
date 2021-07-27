in_bam=$1
genome_fasta=$2
merged_overlap=$3
blacklist=$4
outdir=$5

TOBIAS ATACorrect --bam $in_bam --genome $genome_fasta --peaks $merged_overlap --blacklist $blacklist --outdir $outdir --cores 8 --read_shift 4 -4

