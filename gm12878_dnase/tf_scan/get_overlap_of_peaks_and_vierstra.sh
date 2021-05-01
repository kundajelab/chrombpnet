## 50 bp around summit 
#intersect with motif archetypes
#bedtools intersect -wo -a GM12878.idr.peaks.50bp.around.summit.bed -b /mnt/data/vierstra_motifs/hg38.archetype_motifs.v1.0.bed.gz > GM12878.idr.peaks.50bp.around.summit.overlap.tf.bed

#intersect with individual motifs
#bedtools intersect -nobuf -wo -b GM12878.idr.peaks.50bp.around.summit.bed -a /mnt/data/vierstra_motifs/hg38.all_motifs.v1.0.bed > GM12878.idr.peaks.50bp.around.summit.overlap.tf.all.bed


#intersect with individual motifs, looking at CTCF only
#bedtools intersect  -nobuf -wo -a GM12878.idr.peaks.50bp.around.summit.bed -b /mnt/data/vierstra_motifs/CTCF.HUMAN.hg38.all_motifs.v1.0.bed > CTCF.GM12878.idr.peaks.50bp.around.summit.overlap.tf.all.bed


## NOT LIMITING TO 50 BP AROUND SUMMIT
#bedtools intersect -wo -a GM12878.idr.optimal_peak.narrowPeak.gz -b /mnt/data/vierstra_motifs/hg38.archetype_motifs.v1.0.bed.gz > GM12878.idr.peaks.overlap.tf.archetypes.bed
#split by TSS+/-1kb and other
bedtools intersect -a GM12878.idr.peaks.overlap.tf.archetypes.bed -b gencode.v31.hg38.tss.bed > GM12878.idr.peaks.overlap.tf.archetypes.TSS.1kb.bed
bedtools intersect -v -a GM12878.idr.peaks.overlap.tf.archetypes.bed -b gencode.v31.hg38.tss.bed > GM12878.idr.peaks.overlap.tf.archetypes.NOT.TSS.1kb.bed


#bedtools intersect -wo -b GM12878.idr.optimal_peak.narrowPeak.gz -a /mnt/data/vierstra_motifs/hg38.all_motifs.v1.0.bed.gz > GM12878.idr.peaks.overlap.tf.all.bed
