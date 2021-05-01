#intersect with individual motifs
bedtools intersect -nobuf -wo -b HEPG2.idr.peaks.50bp.around.summit.bed -a /mnt/data/vierstra_motifs/hg38.all_motifs.v1.0.bed > HEPG2.idr.peaks.50bp.around.summit.overlap.tf.all.bed
