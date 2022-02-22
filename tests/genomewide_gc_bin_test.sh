#!/bin/bash

#pull down chr21 of hg38 for test
wget -O chr21.fa.gz https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz 
gzip -df chr21.fa.gz
echo "downloaded chr21 hg38 fasta file" 

#test for several strides
for stride in 100 500 1000
do
    echo $stride 
    python ../src/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py -g chr21.fa -o $stride.chr21.gc.bed -f 2114 -s $stride 
done

#assert we get the expectd md5sum values
#2911f95c85673d86cff8daa4eaca5250  100.chr21.gc.bed
#0bec3a55b97a2fa1a4f62e5abb375653  1000.chr21.gc.bed
#43c3d7185c4748532b36959f3c0ca346  500.chr21.gc.bed

#get md5sum values for the genomewide gc content
md5sum_100=($(md5sum 100.chr21.gc.bed))
md5sum_500=($(md5sum 500.chr21.gc.bed))
md5sum_1000=($(md5sum 1000.chr21.gc.bed))
echo $md5sum_100
echo $md5sum_500
echo $md5sum_1000

if [ $md5sum_100 = "2911f95c85673d86cff8daa4eaca5250" ]; 
then
    echo "obtained expected md5sum for 100.chr21.gc.bed"
else
    echo "md5sum $md5sum_100 does not match expected value 2911f95c85673d86cff8daa4eaca5250"
    exit 1
fi

if [ $md5sum_500 = "43c3d7185c4748532b36959f3c0ca346" ];
then
    echo "obtained expected md5sum for 500.chr21.gc.bed"
else
    echo "md5sum $md5sum_500 does not match expected value 43c3d7185c4748532b36959f3c0ca346"
    exit 1
fi

if [ $md5sum_1000 = "0bec3a55b97a2fa1a4f62e5abb375653" ]; 
then
    echo "obtained expected md5sum for 100.chr21.gc.bed"
else
    echo "md5sum $md5sum_100 does not match expected value 0bec3a55b97a2fa1a4f62e5abb375653"
    exit 1
fi


