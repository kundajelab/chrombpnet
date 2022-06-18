#!/bin/bash
set -e
set -o pipefail
set -u

# Converts a set of bams to + and - strand 5' end bigwigs used for training BPNet
# By default keeps duplicates, change samtools view flag to change.

# make sure to give a prefix, final files are named ${OUTPREFIX}_{plus,minus}.bw
OUTPREFIX=$1  # e.g. /path/to/dir/prefix
REFCHROMSZ=$2 # e.g. /path/to/genome.sizes.txt
INTA=$3 # tagAlign file

if ! [ -x "$(command -v bedtools)" ] ; then
    echo "Bedtools not found" >&2
    exit 1
fi

if ! [ -x "$(command -v bedGraphToBigWig)" ] ; then
    echo "bedGraphToBigWig not found" >&2
    exit 1
fi

printf "Making BedGraphs\n"
# tagAlign -> bedGraph
# genomecov documentation mentions only need to sort by chromosome
# bedtools genomecov -5 -bg -g $REFCHROMSZ -i <(zcat $INTA | awk -v OFS='\t' '{if ($6=="+") {print $0} else {print $1,$2+1,$3+1,$4,$5,$6}}' | sort -k1,1) | bedtools sort  > ${OUTPREFIX}.bedGraph

bedtools genomecov -5 -bg -g $REFCHROMSZ -i <(zcat $INTA | awk -v OFS='\t' '{if ($6=="+") {print $0} else {print $1,$2+1,$3+1,$4,$5,$6}}') > ${OUTPREFIX}-unsorted.bedgraph

printf "Sorting BedGraphs\n"
sort -k1,1 -k2,2n ${OUTPREFIX}-unsorted.bedgraph > ${OUTPREFIX}.bedGraph

printf "Making BigWigs\n"
# bedGraph -> BigWig
bedGraphToBigWig ${OUTPREFIX}.bedGraph $REFCHROMSZ ${OUTPREFIX}.bw

# remove intermediates
rm ${OUTPREFIX}.bedGraph
${OUTPREFIX}-unsorted.bedgraph
