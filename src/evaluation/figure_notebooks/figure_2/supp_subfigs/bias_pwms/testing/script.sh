#!/bin/bash
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

# DNASE is SE data
# ATAC is PE
## download ATAC bam
wget https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam -O atac_rep1.bam
wget https://www.encodeproject.org/files/ENCFF860XAE/@@download/ENCFF860XAE.bam -O dnase_rep1.bam

samtools sort atac_rep1.bam -o atac_rep1.sorted.bam
samtools sort dnase_rep1.bam -o dnase_rep1.sorted.bam

dnase_in_bam=dnase_rep1.sorted.bam
atac_in_bam=atac_rep1.sorted.bam

wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O hg38.fa.gz
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O hg38.chrom.sizes

gunzip hg38.fa.gz
bgzip hg38.fa

genome=hg38.fa.gz
chrom_sizes=hg38.chrom.sizes


# subsample only "chr20" for testing

samtools index $atac_in_bam
samtools index $dnase_in_bam

samtools view -b $atac_in_bam chr20 > out.bam

## no shift ATAC data  
echo "Generate images for atac with no shift on positive strand"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | bedtools genomecov -bg -5 -strand + -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes atac_pos.bw
rm tmp2
python build_pwm_from_bigwig.py -i  atac_pos.bw -g $genome -o atac_no_shift_pos.pdf -c "chr20" -cz $chrom_sizes
rm atac_pos.bw

echo "Generate images for atac with no shift on negative strand"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | bedtools genomecov -bg -5 -strand - -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes atac_neg.bw
rm tmp2
python build_pwm_from_bigwig.py -i  atac_neg.bw -g $genome -o atac_no_shift_neg.pdf -c "chr20" -cz $chrom_sizes
rm atac_neg.bw


## shift ATAC data by +4/-4
echo "Generate images for dnase with +4/-4 shift"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-4,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o atac_44_shift.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw

samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-4,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -strand + -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.pos.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.pos.bw -g $genome -o atac_44_shift_pos.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.pos.bw


samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-4,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -strand - -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.neg.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.neg.bw -g $genome -o atac_44_shift_neg.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.neg.bw


## shift ATAC data by +4/-5
echo "Generate images for dnase with +4/-5 shift"
samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-5,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o atac_45_shift.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw


samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-5,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -strand + -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.pos.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.pos.bw -g $genome -o atac_45_shift_pos.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.pos.bw

samtools view -b -@50 out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3-5,$4,$5,$6} else if ($6=="+") {print $1,$2+4,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -strand - -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.neg.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.neg.bw -g $genome -o atac_45_shift_neg.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.neg.bw

rm out.bam



# subsample only "chr20" for testing

samtools view -b $dnase_in_bam chr20 > out.bam

## no shift DNASE data  
echo "Generate images for dnase with no shift"

samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o dnase_no_shift.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw


samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | bedtools genomecov -bg -5 -i stdin -strand + -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o dnase_no_shift_pos.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw

samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | bedtools genomecov -bg -5 -strand -  -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o dnase_no_shift_neg.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw



## shift DNASE with 0/1
echo "Generate images for dnase with 0/1 shift"
samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o dnase_01_shift.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw


samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -strand +  -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o dnase_01_shift_pos.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw


samtools view -b -@50 -F780 -q30  out.bam | bedtools bamtobed -i stdin | awk -v OFS="\t" '{if ($6=="-"){print $1,$2,$3+1,$4,$5,$6} else if ($6=="+") {print $1,$2,$3,$4,$5,$6}}' | bedtools genomecov -bg -5 -i stdin -strand -  -g $chrom_sizes | bedtools sort -i stdin > tmp2
bedGraphToBigWig tmp2 $chrom_sizes shifted.unstranded.bw
rm tmp2
python build_pwm_from_bigwig.py -i  shifted.unstranded.bw -i shifted.unstranded.bw -g $genome -o dnase_01_shift_neg.pdf -c "chr20" -cz $chrom_sizes 
rm shifted.unstranded.bw


rm out.bam


