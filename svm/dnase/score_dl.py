#calculates auPRC on positive and negative inputs to lsgkm prediction
import argparse
import pandas as pd
import pdb
from sklearn.metrics import average_precision_score, precision_recall_curve


task_name_to_colname={}
task_name_to_colname['dnase_c']=0
task_name_to_colname['dnase_v']=1
task_name_to_colname['sw480']=2
task_name_to_colname['hct116']=3
task_name_to_colname['colo205']=4

colname_to_task_name={}
colname_to_task_name[0]='dnase_c'
colname_to_task_name[1]='dnase_v'
colname_to_task_name[2]='sw480'
colname_to_task_name[3]='hct116'
colname_to_task_name[4]='colo205'

def parse_args():
    parser=argparse.ArgumentParser(description="calculates auPRC on hdf5 prediction files from deep learning models")
    parser.add_argument("--labels_hdf5",nargs="+")
    parser.add_argument("--predictions_hdf5",nargs="+")
    parser.add_argument("--outf")
    parser.add_argument("--multitask",default=False,action="store_true")
    return parser.parse_args()

def score_one_task(preds,labels):
    pass

def main():
    args=parse_args()
    #make sure the label hdf5 inputs are matched with prediction hdf5 inputs
    assert(len(args.labels_hdf5)==len(args.predictions_hdf5))    
    num_datasets=len(args.labels_hdf5)
    sample_to_auprc=dict()
    for i in range(num_datasets):
        print(args.labels_hdf5[i])
        print(args.predictions_hdf5[i])
        cur_preds=pd.read_hdf(args.predictions_hdf5[i])
        cur_labels=pd.read_hdf(args.labels_hdf5[i])
        num_tasks=cur_preds.shape[1]
        if (num_tasks > 1) and (args.multitask==True):
            #score all the tasks
            for cur_task in range(num_tasks):
                task_labels=cur_labels[cur_task]
                task_preds=cur_preds[cur_task]
                cur_subset=pd.concate([task_labels,task_preds],axis=1).dropna()
                cur_subset.columns=['labels','preds']
                task_name=colname_to_task_name[cur_task]
                cur_sample=args.labels_hdf5[i].strip('.labels.0')+'.'+task_name
                cur_auprc=average_precision_score(cur_subset['labels'],cur_subset['preds'])
                sample_to_auprc[cur_sample]=cur_auprc
        elif (num_tasks > 1):
            #get the actual task column
            for key in task_name_to_colname:
                if key in args.labels_hdf5[i]:
                    #extract the corresponding column
                    cur_task_colname=task_name_to_colname[key]
                    cur_labels=cur_labels[cur_task_colname]
                    #assert the labels and predictions dataframes are matched
                    assert key in args.predictions_hdf5[i] 
                    cur_preds=cur_preds[cur_task_colname]
                        
                    cur_data=pd.concat((cur_preds,cur_labels),axis=1)
                    cur_data=cur_data.dropna()
                    cur_data.columns=['preds','labels']
                    cur_auprc=average_precision_score(cur_data['labels'],cur_data['preds'])
                    cur_sample=args.labels_hdf5[i].strip('.labels.0') 
                    sample_to_auprc[cur_sample]=cur_auprc
        else:
            cur_data=pd.concat((cur_preds,cur_labels),axis=1).dropna()
            cur_data.columns=['preds','labels']
            #pdb.set_trace() 
            cur_auprc=average_precision_score(cur_data['labels'],cur_data['preds'])
            cur_sample=args.labels_hdf5[i].strip('.labels.0') 
            sample_to_auprc[cur_sample]=cur_auprc
    print(sample_to_auprc)
    outf=open(args.outf+"/perf.metrics.txt",'w')
    outf.write('Dataset\tauPRC\n')
    for key in sample_to_auprc:
        outf.write(key+'\t'+str(sample_to_auprc[key])+'\n')
    outf.close()
    
if __name__=="__main__":
    main()
    
