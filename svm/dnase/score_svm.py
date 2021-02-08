#calculates auPRC on positive and negative inputs to lsgkm prediction
import argparse
import pandas as pd
import pdb
from sklearn.metrics import average_precision_score, precision_recall_curve

def parse_args():
    parser=argparse.ArgumentParser(description="calculates auPRC on positive and negative inputs to lsgkm prediction")
    parser.add_argument('--predictions_on_pos_regions',nargs="+")
    parser.add_argument('--predictions_on_neg_regions',nargs="+")
    parser.add_argument('--outf')
    return parser.parse_args()

def main():
    args=parse_args()
    #outf=open(args.outf,'w')
    num_datasets=len(args.predictions_on_pos_regions)
    sample_to_auprc=dict() 
    for i in range(num_datasets):
        print(args.predictions_on_pos_regions[i])
        print(args.predictions_on_neg_regions[i])
        pos_region_preds=pd.read_csv(args.predictions_on_pos_regions[i],header=None,sep='\t')
        neg_region_preds=pd.read_csv(args.predictions_on_neg_regions[i],header=None,sep='\t')
        pos_region_preds['label']=1
        pos_region_preds.columns=['region','pred','label']
        neg_region_preds['label']=-1
        neg_region_preds.columns=['region','pred','label']
        #concatenate positive & negative test sets
        all_preds= pd.concat([pos_region_preds, neg_region_preds],axis=0).dropna()
        #print(all_preds.head())
        #calculate the auprc
        cur_auprc=average_precision_score(all_preds['label'],all_preds['pred'])
        cur_sample=args.predictions_on_pos_regions[i].strip('.positives') 
        sample_to_auprc[cur_sample]=cur_auprc
        #get values for precision recall curve
    print(sample_to_auprc)
    outf=open(args.outf+"perf.metrics.txt",'w')
    outf.write('Dataset\tauPRC\n')
    for key in sample_to_auprc:
        outf.write(key+'\t'+str(sample_to_auprc[key])+'\n')
    outf.close()
    
    

    

if __name__=="__main__":
    main()
    
