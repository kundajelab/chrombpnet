peak_bed=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/nautilus_runs_jun16/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/interpret/merged.HEPG2.interpreted_regions.bed
out_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2//HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/moods/
in_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2//HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/
fasta=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
mode=counts
mkdir $out_dir/$mode
echo $out_dir

FILE=$out_dir/merged_peaks_no_lacklist.w1000.bed
if [ -f "$FILE" ]; then
    echo "$FILE exists."
else 
    echo "$FILE does not exist."

    bedtools slop -i /mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz -g /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes -b 1057 > temp.bed
    bedtools intersect -v -a $peak_bed -b temp.bed | awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' > $out_dir/merged_peaks_no_lacklist.w1000.bed
    bed_file=$out_dir/merged_peaks_no_lacklist.w1000.bed
    bedtools getfasta -fi $fasta -bed $bed_file -fo $out_dir/merged.peaks.fasta
    rm temp.bed
fi

va=$(grep -v "^>" $out_dir/merged.peaks.fasta | cat |  python3 -c "import sys; print(sum(line.upper().count('A') for line in sys.stdin))")
vc=$(grep -v "^>" $out_dir/merged.peaks.fasta | cat |  python3 -c "import sys; print(sum(line.upper().count('C') for line in sys.stdin))")
vg=$(grep -v "^>" $out_dir/merged.peaks.fasta | cat |  python3 -c "import sys; print(sum(line.upper().count('G') for line in sys.stdin))")
vt=$(grep -v "^>" $out_dir/merged.peaks.fasta | cat |  python3 -c "import sys; print(sum(line.upper().count('T') for line in sys.stdin))")

tot=$(($va + $vc + $vg + $vt))
echo $tot
vaf=$(python -c "print($va / $tot)")
vgf=$(python -c "print($vg / $tot)")
vcf=$(python -c "print($vc / $tot)")
vtf=$(python -c "print($vg / $tot)")

echo $vaf
echo $vcf
echo $vgf
echo $vtf


for pfm in $in_dir/ppms/$mode/*
    do
        base=$(basename -- "$pfm")
        tf="${base%.*}"
        tf="${tf%.*}"
        tf="${tf%.*}"
        echo $tf
        moods-dna.py -m $pfm -s $out_dir/merged.peaks.fasta -p 0.0001 --lo-bg $vaf $vcf $vgf $vtf --batch > $out_dir/$mode/$tf.hits.csv &
   done

wait

cat $out_dir/$mode/metacluster*.hits.csv > $out_dir/$mode/all_metaclusters.all_patterns.hits.csv




