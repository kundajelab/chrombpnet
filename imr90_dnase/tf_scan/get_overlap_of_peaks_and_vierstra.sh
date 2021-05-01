#intersect with individual motifs
bedtools intersect -nobuf -wo -b IMR90.idr.peaks.50bp.around.summit.bed -a /mnt/data/vierstra_motifs/hg38.all_motifs.v1.0.bed > IMR90.idr.peaks.50bp.around.summit.overlap.tf.all.bed
