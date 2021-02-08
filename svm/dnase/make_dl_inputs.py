#generate inputs for predictions by gc-corrected dl multi-tasked models
import argparse
import pandas as pd
import numpy as np 
import pdb

def parse_args():
    parser=argparse.ArgumentParser(description="form test inputs for dl")
    parser.add_argument("--tasks",nargs="+")
    parser.add_argument("--positives",nargs="+")
    parser.add_argument("--negatives",nargs="+")
    parser.add_argument("--store_gc",action="store_true",default=False)
    parser.add_argument('--out_prefix') 
    return parser.parse_args()

def main():
    args=parse_args()
    num_tasks=len(args.tasks)
    df_dict=dict() 
    for i in range(len(args.tasks)):
        cur_task=args.tasks[i]
        cur_pos=pd.read_csv(args.positives[i],header=None,sep='\t')
        cur_neg=pd.read_csv(args.negatives[i],header=None,sep='\t')
        cur_pos[cur_task]=1
        cur_neg[cur_task]=0
        labels=cur_pos[[0,1,2,cur_task]].append(cur_neg[[0,1,2,cur_task]])
        for task in args.tasks:
            if task!=cur_task:
                labels[task]=np.nan
        labels['CHR']=labels[0]
        labels['START']=labels[1]
        labels['END']=labels[2]
        labels.set_index(['CHR','START','END'])
        labels.drop(columns=[0,1,2],inplace=True)
        col_order=['CHR','START','END']+args.tasks
        labels=labels[col_order]
        print(labels.head())
        labels.to_hdf(args.out_prefix+'.'+cur_task+'.labels.hdf5',key='data',mode='w',format='table',min_itemsize={'CHR':30})
        gc=cur_pos[[0,1,2,3]].append(cur_neg[[0,1,2,3]])
        gc.columns=['CHR','START','END','gc']
        gc.set_index(['CHR','START','END'])
        print(gc.head())
        gc.to_hdf(args.out_prefix+'.'+cur_task+'.gc.hdf5',key='data',mode='w',format='table',min_itemsize={'CHR':30})
        

        

if __name__=="__main__":
    main()
    
    
