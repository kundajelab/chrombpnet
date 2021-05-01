## NOT LIMITING TO 50 BP AROUND SUMMIT
#bedtools intersect -wo -a K562.dnase.idr.optimal_peak.narrowPeak.gz -b /mnt/data/vierstra_motifs/hg38.all_motifs.v1.0.bed.gz > K562.idr.peaks.overlap.tf.all.bed
#bedtools intersect -wo -a K562.dnase.idr.optimal_peak.narrowPeak.gz -b /mnt/data/vierstra_motifs/hg38.archetype_motifs.v1.0.bed.gz > K562.idr.peaks.overlap.tf.archetypes.bed
for chrom in X Y #`seq 1 22` X Y
do
    bedtools intersect -wo -a K562.dnase.idr.optimal_peak.narrowPeak.gz -b /mnt/data/vierstra_motifs/chr$chrom.hg38.all_motifs.v1.0.bed.gz > chr$chrom.K562.idr.peaks.overlap.tf.all.bed &
    echo $chrom 
done
