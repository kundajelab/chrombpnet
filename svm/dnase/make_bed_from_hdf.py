#convert hdf5 label files into bed format 
import sys
task=sys.argv[1]
fold=sys.argv[2]
import pandas as pd
in_prefix="/srv/scratch/annashch/5_cell_lines_bias_correction/gc_covariate/classification/"
out_prefix="/srv/scratch/annashch/5_cell_lines_bias_correction/svm"
data=pd.read_hdf(in_prefix+"/"+task+"/"+"predictions.DNASE."+task+".classificationlabels.withgc."+fold+".labels.0")
data.to_csv(out_prefix+"/"+task+"/"+"svm_predictions_svmtrainset_genometestset"+"/"+"labels."+str(fold)+".bed",sep='\t',index=True,header=False)
