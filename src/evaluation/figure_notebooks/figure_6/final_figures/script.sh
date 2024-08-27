#cp /mnt/lab_data3/anusri/histone_expts/all_qtl_analysis/gm12878_sequence_sets/test_set/deltasv,/41588_2015_BFng3331_MOESM26_ESM.csv dsqtls/
#cp /mnt/lab_data2/anusri/chrombpnet/results/variant_data/dsqtls/process/fetch_enformer/enformer_predictions.tsv dsqtls/
#zcat /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/GSE31388_dsQtlTable.txt.gz > GSE31388_dsQtlTable.txt

cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/ATAC/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR637XSC.ATAC.GM12878.fold.mean.scores.tsv
dd=250M
cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/ATAC_$dd/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR637XSC.$dd.ATAC.GM12878.fold.mean.scores.tsv
dd=100M
cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/ATAC_$dd/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR637XSC.$dd.ATAC.GM12878.fold.mean.scores.tsv
dd=50M
cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/ATAC_$dd/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR637XSC.$dd.ATAC.GM12878.fold.mean.scores.tsv
dd=25M
cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/ATAC_$dd/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR637XSC.$dd.ATAC.GM12878.fold.mean.scores.tsv
dd=5M
cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/ATAC_$dd/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR637XSC.$dd.ATAC.GM12878.fold.mean.scores.tsv
cp /mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/DNASE/summary.mean.variant_scores_new_2.tsv dsqtls/ENCSR000EMT.DNASE.GM12878.fold.mean.scores.tsv




