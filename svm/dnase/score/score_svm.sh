#default kernel
#for fold in `seq 0 9`
#do
#    for task in C_merged V_merged sw480 hct116 colo205
#    do
#	kerasAC_score --labels_hdf5 SVM_predictions_labels_df/$task.$fold.labels.hdf5 \
#		      --predictions_hdf5 SVM_predictions_labels_df/$task.$fold.predictions.hdf5 \
#		      --performance_metrics_classification_file SVM_perform/perf.$task.$fold \
#		      --task_indices 0
#    done
#done

#rbf kernel
for fold in `seq 0 9`
do
    for task in C_merged V_merged sw480 hct116 colo205
    do
	kerasAC_score --labels_hdf5 SVM_predictions_labels_df_rbf/$task.$fold.labels.hdf5 \
		      --predictions_hdf5 SVM_predictions_labels_df_rbf/$task.$fold.predictions.hdf5 \
		      --performance_metrics_classification_file SVM_perform_rbf/perf.$task.$fold \
		      --task_indices 0
    done
done

