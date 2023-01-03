bed=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/data/30K.subsample.overlap.bed 


bedtools sort -i $bed | bedtools merge -i stdin > merged.bed
bedtools sort -i $bed | bedtools merge -i stdin | bedtools intersect -a /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/k562_tf_chip_combine.bed -b stdin -wa -f 1.0 | bedtools sort -i stdin | uniq > chip_in_bed.bed

# not needed below step
#bedtools intersect -a modisco_hits_k562_atac_uncorrected.bed -b chip_in_bed.bed -wa -f 0.5 > modisco_chip.intersect

