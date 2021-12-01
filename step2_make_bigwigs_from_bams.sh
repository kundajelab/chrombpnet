in_bam=$1
bigwig_dir=$2
data_type=$3
reference_fasta=$4
chrom_sizes=$5

bash $PWD/src/helpers/preprocessing/bam_to_bigwig.sh $in_bam $bigwig_dir $data_type $chrom_sizes 
python $PWD/src/helpers/preprocessing/analysis/build_pwm_from_bigwig.py -i $bigwig_dir/shifted.sorted.bam.chrombpnet.unstranded.bw  -g $reference_fasta -o $bigwig_dir/bias_pwm.png -c "chr20" -cz $chrom_sizes 
