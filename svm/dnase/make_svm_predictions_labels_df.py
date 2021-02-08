#transform lsgkm outputs to labels/prediction dataframes for scoring with kerasAC
import argparse
import pandas as pd
import pdb

def parse_args():
    parser=argparse.ArgumentParser(description="Convert SVM predictions to a data frame")
    parser.add_argument("--pos_predictions")
    parser.add_argument("--neg_predictions")
    parser.add_argument("--out_prefix")
    return parser.parse_args()

def main():
    args=parse_args()
    pos_predictions=pd.read_csv(args.pos_predictions,header=None,sep='\t')
    #pos_predictions=pos_predictions[0].str.split('_',expand=True)
    neg_predictions=pd.read_csv(args.neg_predictions,header=None,sep='\t')
    #neg_predictions=neg_predictions[0].str.split('_',expand=True)
    pos_labels=pos_predictions.copy()
    neg_labels=neg_predictions.copy()
    pos_labels[1]=1
    neg_labels[1]=0
    #concatenate
    labels=pd.concat([pos_labels,neg_labels],axis=0)
    predictions=pd.concat([pos_predictions,neg_predictions],axis=0)
    predictions[1]=1.0*pd.to_numeric(predictions[1]>0)
    labels.columns=['CHR_START_END',0]
    predictions.columns=['CHR_START_END',0]
    labels=labels.set_index(['CHR_START_END'])
    predictions=predictions.set_index(['CHR_START_END'])
    #labels[0]=pd.to_numeric(labels[0],errors='coerce')
    #pdb.set_trace() 
    #predictions[0]=(pd.to_numeric(predictions[0],errors='coerce')>0)*1.0
    
    
    #save to hdf5
    labels.to_hdf(args.out_prefix+'.labels.hdf5',key='data',mode='w',format='table')
    predictions.to_hdf(args.out_prefix+'.predictions.hdf5',key='data',mode='w',format='table')

if __name__=="__main__":
    main()
    
