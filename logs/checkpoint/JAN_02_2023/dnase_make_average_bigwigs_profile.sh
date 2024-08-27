oakdir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"
celltype="HEPG2"
zcat  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.bed" >  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
wc -l  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
python  /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-h5 $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.h5" \
	-r  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed" \
	-c /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes \
	-o  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.bw" \
	-s  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.stat" \
	-t 1

celltype="K562"
zcat  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.bed" >  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
wc -l  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
python  /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-h5 $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.h5" \
	-r  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed" \
	-c /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes \
	-o  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.bw" \
	-s  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.stat" \
	-t 1

celltype="GM12878"
zcat  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.bed" >  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
wc -l  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
python  /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-h5 $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.h5" \
	-r  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed" \
	-c /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes \
	-o  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.bw" \
	-s  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.stat" \
	-t 1


celltype="IMR90"
zcat  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.bed" >  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
wc -l  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
python  /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-h5 $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.h5" \
	-r  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed" \
	-c /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes \
	-o  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.bw" \
	-s  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.stat" \
	-t 1


celltype="H1ESC"
zcat  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.bed" >  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
wc -l  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed"
python  /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-h5 $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.h5" \
	-r  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores_new_compressed.unzip.bed" \
	-c /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes \
	-o  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.bw" \
	-s  $oakdir/$celltype/merge_folds_new/$celltype"_folds_merged.profile_scores.stat" \
	-t 1
