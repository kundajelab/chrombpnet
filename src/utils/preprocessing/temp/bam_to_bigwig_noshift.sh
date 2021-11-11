#!/bin/bash

input_bam=$1
output_dir=$2
data_type=$3
chrom_sizes=$4

if [ "$data_type" = "DNASE_SE" ] ; then
    echo "shift DNASE SE data"
    samtools view -b -@50 -F780 -q30  $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_dir/tmp2
    bedGraphToBigWig $output_dir/tmp2 $chrom_sizes $output_dir/shifted.sorted.bam.chrombpnet.unstranded.bw
    rm $output_dir/tmp2
elif [ "$data_type" = "DNASE_PE" ] ; then
    echo "shift DNASE PE data"
    samtools view -b -@50 $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_dir/tmp2
    bedGraphToBigWig $output_dir/tmp2 $chrom_sizes $output_dir/shifted.sorted.bam.chrombpnet.unstranded.bw
    rm $output_dir/tmp2
elif [ "$data_type" = "ATAC_PE" ] ; then
    echo "shift ATAC PE data"
    samtools view -b -@50 $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_dir/tmp2
    bedGraphToBigWig $output_dir/tmp2 $chrom_sizes $output_dir/shifted.sorted.bam.chrombpnet.unstranded.bw
    rm $output_dir/tmp2
elif [ "$data_type" = "ATAC_SE" ] ; then
    echo "shift ATAC SE data"
    samtools view -b -@50 -F780 -q30  $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_dir/tmp2
    bedGraphToBigWig $output_dir/tmp2 $chrom_sizes $output_dir/shifted.sorted.bam.chrombpnet.unstranded.bw
    rm $output_dir/tmp2
else
    echo "unknown data type"
    echo $data_type
fi
