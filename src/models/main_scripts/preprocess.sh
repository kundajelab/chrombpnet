#!/bin/bash

## WARNING - Check samtools flags when creating bigwigs
## WARNING - Check tn5 shift

in_bam=$1
interm=$2
samtools_flag=$3
is_filtered=$4
type=$5
neg_shift=$6
chrom_sizes=$7

echo $chrom_sizes
## shift bam files

if [ "$type" = DNASE ] ; then
    echo "shift DNASE data"
    samtools view -b -@50 -F780 -q30  $in_bam | bedtools bamtobed -i stdin > $interm/data.bed
    awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2} else if ($6 == "-") {$3 = $3 + 1} print $0}' $interm/data.bed > $interm/shifted.bed
    bedToBam -i $interm/shifted.bed -g  $chrom_sizes  > $interm/shifted_0_1.bam
    shifted_bam=$interm/shifted_0_1
else
    echo "shift ATAC data"
    alignmentSieve -p 56 -b $in_bam -o $interm/shifted_4_$neg_shift.bam --shift 4 -$neg_shift $neg_shift -4  --filterMetrics $interm/log_4_$neg_shift.txt
    shifted_bam=$interm/shifted_4_$neg_shift
fi

# get unstranded bigwigs

samtools sort -@ 16 $shifted_bam.bam -o $shifted_bam.sorted.bam


cur_bam=$shifted_bam.sorted.bam
samtools index $cur_bam

if [ "$is_filtered" = True ] ; then
    echo "not filtering data"
    bamCoverage -p16 -v --binSize 1  --Offset 1 1 -b $cur_bam -o $interm/shifted.sorted.bam.bpnet.unstranded.bw
else
    echo "filtering data"
    bamCoverage -p16 -v --binSize 1 --samFlagExclude $samtools_flag  --Offset 1 1 --minMappingQuality 30 -b $cur_bam -o $interm/shifted.sorted.bam.bpnet.unstranded.bw
fi
