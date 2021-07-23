#!/bin/bash

## WARNING - Check samtools flags when creating bigwigs
## WARNING - Check tn5 shift

in_bam=$1
interm=$2
samtools_flag=$3
is_filtered=$4
type=$5

## shift bam files

if [ "$type" = DNASE ] ; then
    echo "shift DNASE data"
    alignmentSieve -p 56 -b $in_bam -o $interm/shifted_0_1.bam --shift 0 1 --filterMetrics $interm/log_0_1.txt
    shifted_bam=$interm/shifted_0_1
else
    echo "shift ATAC data"
    alignmentSieve -p 56 -b $in_bam -o $interm/shifted_4_4.bam --shift 4 -4 4 -4  --filterMetrics $interm/log_4_4.txt
    shifted_bam=$interm/shifted_4_4
fi

shifted_bam=$interm/shifted_0_1

# get unstranded bigwigs

samtools sort -@ 16 $shifted_bam.bam -o $shifted_bam.sorted.bam


cur_bam=$shifted_bam.sorted.bam
samtools index $cur_bam

if [ "$is_filtered" = True ] ; then
    echo "not filtering data"
    bamCoverage -p16 -v --binSize 1  --Offset 1 1 --minMappingQuality 30 -b $cur_bam -o $interm/shifted_4_4.sorted.bam.bpnet.unstranded.bw
else
    echo "filtering data"
    bamCoverage -p16 -v --binSize 1 --samFlagExclude $samtools_flag  --Offset 1 1 --minMappingQuality 30 -b $cur_bam -o $interm/shifted_4_4.sorted.bam.bpnet.unstranded.bw
fi
