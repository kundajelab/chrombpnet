#!/bin/bash
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

# DNASE is SE data
# ATAC is PE
dnase_in_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/unfiltered_bams/GM12878/GM12878.unfiltered.sorted.bam
#/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam
atac_in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes


# subsample only "chr20" for testing

samtools view -b $dnase_in_bam chr20 > out.bam

## no shift DNASE data  

echo "Generate images for dnase with no shift on plus strand"
samtools view -b -@50 -F796 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes unshifted.sorted.bam.chrombpnet.plus.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  unshifted.sorted.bam.chrombpnet.plus.bw -i unshifted.sorted.bam.chrombpnet.plus.bw -g $genome -o dnase_no_shift_plus.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 

echo "Generate images for dnase with no shift on negative strand"
samtools view -b -@50 -F780 -f16  -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes unshifted.sorted.bam.chrombpnet.minus.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  unshifted.sorted.bam.chrombpnet.minus.bw -i unshifted.sorted.bam.chrombpnet.minus.bw -g $genome -o dnase_no_shift_minus.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 


echo "Generate images for dnase with no shift"
samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes unshifted.sorted.bam.chrombpnet.unstranded.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  unshifted.sorted.bam.chrombpnet.unstranded.bw -i unshifted.sorted.bam.chrombpnet.unstranded.bw -g $genome -o dnase_no_shift.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 


## shift DNASE with 0/1
echo "Generate images for dnase with no shift on plus strand"
samtools view -b -@50 -F796 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.sorted.bam.chrombpnet.plus.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  shifted.sorted.bam.chrombpnet.plus.bw -i shifted.sorted.bam.chrombpnet.plus.bw -g $genome -o dnase_01_shift_plus.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 

echo "Generate images for dnase with no shift on negative strand"
samtools view -b -@50 -F780 -f16  -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.sorted.bam.chrombpnet.minus.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  shifted.sorted.bam.chrombpnet.minus.bw -i shifted.sorted.bam.chrombpnet.minus.bw -g $genome -o dnase_01_shift_minus.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 


echo "Generate images for dnase with 0/1 shift"
samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.sorted.bam.chrombpnet.unstranded.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  shifted.sorted.bam.chrombpnet.unstranded.bw -i shifted.sorted.bam.chrombpnet.unstranded.bw -g $genome -o dnase_01_shift.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 

rm out.bam


# subsample only "chr20" for testing

samtools view -b $atac_in_bam chr20 > out.bam

## no shift ATAC data  
echo "Generate images for atac with no shift"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.sorted.bam.chrombpnet.unstranded.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  shifted.sorted.bam.chrombpnet.unstranded.bw -i shifted.sorted.bam.chrombpnet.unstranded.bw -g $genome -o atac_no_shift.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 
rm shifted.sorted.bam.chrombpnet.unstranded.bw


## shift ATAC data by +4/-4
echo "Generate images for dnase with +4/-4 shift"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-4,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.sorted.bam.chrombpnet.unstranded.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  shifted.sorted.bam.chrombpnet.unstranded.bw -i shifted.sorted.bam.chrombpnet.unstranded.bw -g $genome -o atac_44_shift.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 
rm shifted.sorted.bam.chrombpnet.unstranded.bw


## shift ATAC data by +4/-5
echo "Generate images for dnase with +4/-5 shift"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-5,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.sorted.bam.chrombpnet.unstranded.bw
rm tmp2
python ../../analysis/build_pwm_from_bigwig.py -i  shifted.sorted.bam.chrombpnet.unstranded.bw -i shifted.sorted.bam.chrombpnet.unstranded.bw -g $genome -o atac_45_shift.png -c "chr20" -cz /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes 
rm shifted.sorted.bam.chrombpnet.unstranded.bw

rm out.bam
