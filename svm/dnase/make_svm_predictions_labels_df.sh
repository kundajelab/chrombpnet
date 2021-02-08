#default kernel
#for fold in `seq 0 9`
#do
#    for task in C_merged V_merged sw480 hct116 colo205
#    do
#	python make_svm_predictions_labels_df.py \
#	       --pos_predictions SVM_predictions/predictions.$task.$fold.positives \
#	       --neg_predictions SVM_predictions/predictions.$task.$fold.negatives \
#	       --out_prefix SVM_predictions_labels_df/$task.$fold
#    done
#done

#rbf kernel 
for fold in `seq 0 9`
do
    for task in C_merged V_merged sw480 colo205
    do
	python make_svm_predictions_labels_df.py \
	       --pos_predictions SVM_predictions_rbf/predictions.$task.$fold.positives.gkm.rbf \
	       --neg_predictions SVM_predictions_rbf/predictions.$task.$fold.negatives.gkm.rbf \
	       --out_prefix SVM_predictions_labels_df_rbf/$task.$fold
    done
done
for fold in `seq 0 7` 9
do
    for task in hct116
    do
	python make_svm_predictions_labels_df.py \
	       --pos_predictions SVM_predictions_rbf/predictions.$task.$fold.positives.gkm.rbf \
	       --neg_predictions SVM_predictions_rbf/predictions.$task.$fold.negatives.gkm.rbf \
	       --out_prefix SVM_predictions_labels_df_rbf/$task.$fold
    done
done
