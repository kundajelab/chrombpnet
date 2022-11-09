#!/bin/bash
# exit when any command fails
set -e
set -o pipefail

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'ret=$?;cmd=$last_command;if [ $ret -ne 0 ]; then echo "$cmd failed with error code $ret"; fi' EXIT

input_bam=${1?param missing - input_bam}
output_prefix=${2?param missing - output_prefix}
data_type=${3?param missing - data_type} 
chrom_sizes=${4?param missing - chrom_sizes}
logfile=${5}

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# create the log file
if [ -z "$logfile" ]
  then
    echo "No logfile supplied - creating one"
    logfile=$output_prefix"_bam_to_bigwig.log"
    touch $logfile
fi


if [ "$data_type" = "DNASE_SE" ] ; then
    echo $( timestamp ): "shift DNASE SE data" | tee -a $logfile
    echo $( timestamp ): "samtools view -b -@50 -F780 -q30  $input_bam | bedtools bamtobed -i stdin | awk -v OFS=\"\t\" \'{if (\$6==\"-\"){print \$1,\$2,\$3+1,\$4,\$5,\$6} else if (\$6==\"+\") {print \$1,\$2,\$3,\$4,\$5,\$6}}\' | bedtools genomecov -bg -5 -i stdin -g  $chrom_sizes | bedtools sort -i stdin > $output_prefix.tmp2" | tee -a $logfile
    samtools view -b -@50 -F780 -q30  $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_prefix".tmp2"
    echo $( timestamp ): "bedGraphToBigWig $output_prefix.tmp2 $chrom_sizes $output_prefix\_unstranded.bw" | tee -a $logfile
    bedGraphToBigWig $output_prefix".tmp2" $chrom_sizes $output_prefix"_unstranded.bw"
    echo $( timestamp ): "rm $output_prefix.tmp2" | tee -a $logfile
    rm $output_prefix".tmp2"
elif [ "$data_type" = "DNASE_PE" ] ; then
    echo $( timestamp ): "shift DNASE PE data" | tee -a $logfile
    echo $( timestamp ): "samtools view -b -@50 $input_bam | bedtools bamtobed -i stdin | awk -v OFS=\"\t\" \'{if (\$6==\"-\"){print \$1,\$2,\$3+1,\$4,\$5,\$6} else if (\$6==\"+\") {print \$1,\$2,\$3,\$4,\$5,\$6}}\' | bedtools genomecov -bg -5 -i stdin -g  $chrom_sizes | bedtools sort -i stdin > $output_prefix.tmp2" | tee -a $logfile
    samtools view -b -@50 $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_prefix".tmp2"
    echo $( timestamp ): "bedGraphToBigWig $output_prefix.tmp2 $chrom_sizes $output_prefix\_unstranded.bw" | tee -a $logfile
    bedGraphToBigWig $output_prefix".tmp2" $chrom_sizes $output_prefix"_unstranded.bw"
    echo $( timestamp ): "rm $output_prefix.tmp2" | tee -a $logfile
    rm $output_prefix".tmp2"
elif [ "$data_type" = "ATAC_PE" ] ; then
    echo $( timestamp ): "shift ATAC PE data" | tee -a $logfile
    echo $( timestamp ): "samtools view -b -@50 $input_bam | bedtools bamtobed -i stdin | awk -v OFS=\"\t\" \'{if (\$6==\"-\"){print \$1,\$2,\$3-4,\$4,\$5,\$6} else if (\$6==\"+\") {print \$1,\$2+4,\$3,\$4,\$5,\$6}}\' | bedtools genomecov -bg -5 -i stdin -g  $chrom_sizes | bedtools sort -i stdin > $output_prefix.tmp2" | tee -a $logfile
    samtools view -b -@50 $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-4,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_prefix".tmp2"
    echo $( timestamp ): "bedGraphToBigWig $output_prefix.tmp2 $chrom_sizes $output_prefix\_unstranded.bw" | tee -a $logfile
    bedGraphToBigWig $output_prefix".tmp2" $chrom_sizes $output_prefix"_unstranded.bw"
    echo $( timestamp ): "rm $output_prefix.tmp2" | tee -a $logfile
    rm $output_prefix".tmp2"
elif [ "$data_type" = "ATAC_SE" ] ; then
    echo $( timestamp ): "shift ATAC SE data" | tee -a $logfile
    echo $( timestamp ): "samtools view -b -@50 -F780 -q30  $input_bam | bedtools bamtobed -i stdin | awk -v OFS=\"\t\" \'{if (\$6==\"-\"){print \$1,\$2,\$3-4,\$4,\$5,\$6} else if (\$6==\"+\") {print \$1,\$2+4,\$3,\$4,\$5,\$6}}\' | bedtools genomecov -bg -5 -i stdin -g  $chrom_sizes | bedtools sort -i stdin > $output_prefix.tmp2" | tee -a $logfile
    samtools view -b -@50 -F780 -q30  $input_bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-4,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > $output_prefix".tmp2"
    echo $( timestamp ): "bedGraphToBigWig $output_prefix.tmp2 $chrom_sizes $output_prefix\_unstranded.bw" | tee -a $logfile
    bedGraphToBigWig $output_prefix".tmp2" $chrom_sizes $output_prefix"_unstranded.bw"
    echo $( timestamp ): "rm $output_prefix.tmp2" | tee -a $logfile
    rm $output_prefix".tmp2"
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi
