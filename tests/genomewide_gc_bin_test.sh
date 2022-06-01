#!/bin/bash

#pull down chr21 of hg38 for test
wget -O chr21.fa.gz https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz 
gzip -df chr21.fa.gz
echo "downloaded chr21 hg38 fasta file" 

#test for several strides
for stride in 100 500 1000
do
    echo $stride 
    python ../src/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py -g chr21.fa -o $stride.chr21.gc -f 2114 -s $stride 
done

#assert we get the expectd md5sum values
expected_md5sum_100="1c49d5ff0f5b8e8c09d38360fc50561c"
expected_md5sum_500="f99b40005b496a8f4f7fa9ca04e89bb4"
expected_md5sum_1000="b1c5b8157ff49a9365f0dd6c628afc98"

#get md5sum values for the genomewide gc content
md5sum_100=($(md5sum 100.chr21.gc.bed))
md5sum_500=($(md5sum 500.chr21.gc.bed))
md5sum_1000=($(md5sum 1000.chr21.gc.bed))
echo $md5sum_100
echo $md5sum_500
echo $md5sum_1000

if [ $md5sum_100 = $expected_md5sum_100 ]; 
then
    echo "obtained expected md5sum for 100.chr21.gc.bed"
else
    echo "md5sum $md5sum_100 does not match expected value $expected_md5sum_100"
    exit 1
fi

if [ $md5sum_500 = $expected_md5sum_500 ];
then
    echo "obtained expected md5sum for 500.chr21.gc.bed"
else
    echo "md5sum $md5sum_500 does not match expected value $expected_md5sum_500"
    exit 1
fi

if [ $md5sum_1000 = $expected_md5sum_1000 ]; 
then
    echo "obtained expected md5sum for 1000.chr21.gc.bed"
else
    echo "md5sum $md5sum_100 does not match expected value $expected_md5sum_1000"
    exit 1
fi


