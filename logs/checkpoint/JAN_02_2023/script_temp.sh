oakdir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"
celltype="GM12878"
python  /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
        -h5 $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.counts_scores_new_compressed.h5" \
        -r  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.counts_scores_new_compressed.unzip.bed" \
        -c /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes \
        -o  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.counts_scores.bw" \
        -s  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.counts_scores.stat" \
        -t 1
