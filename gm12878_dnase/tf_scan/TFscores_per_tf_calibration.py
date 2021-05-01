import matplotlib 
from matplotlib import pyplot as plt
import pandas as pd 
import pickle
import pybedtools 
import numpy as np 
from scipy.stats import gmean
import argparse
import pdb

def parse_args():
    parser=argparse.ArgumentParser(description="assign TF deepSHAP scores")
    parser.add_argument("--deepSHAP")
    parser.add_argument("--tf_intersection")
    parser.add_argument("--outf")
    parser.add_argument("--summit_offset",type=int,default=672)
    parser.add_argument("--chrom_col",type=int)
    parser.add_argument("--summit_col",type=int)
    parser.add_argument("--peak_start_col",type=int)
    parser.add_argument("--tf_start_col",type=int)
    parser.add_argument("--tf_end_col",type=int)
    parser.add_argument("--tf_name_col",type=int)
    parser.add_argument("--chrom")
    parser.add_argument("--vierstra_score_col",type=int)
    return parser.parse_args()


def main():
    args=parse_args()
    #load deepshap scores
    deepshap_and_seq=pickle.load(open(args.deepSHAP,'rb'))
    print("loaded deepshap pickle") 
    foreground_shap=deepshap_and_seq['profile_shap']
    foreground_seq=deepshap_and_seq['seq']
    foreground_shap_scaled={} 
    for key in foreground_shap:
        if key[0]==args.chrom:
            foreground_shap_scaled[key]=np.sum(foreground_shap[key]*foreground_seq[key],axis=1)
    print("generated deepSHAP dict for chrom:"+str(args.chrom))

    #load peak intersections with TF    
    intersections=pd.read_csv(args.tf_intersection,header=None,sep='\t',error_bad_lines=False)
    print("loaded intersection file") 
    outf=open(args.outf,'w')
    outf.write('\t'.join(['Chrom','TFStart','TFEnd','TFName','VierstraScore','MeanDeepSHAP','MedianDeepSHAP','90PercentileDeepSHAP','MaxDeepSHAP','MinDeepSHAP'])+'\n')
    count=0
    num_rows=str(intersections.shape[0])
    for index,row in intersections.iterrows():
        count+=1
        if count%10000==0:
            print(str(count)+'/'+num_rows)
        chrom=row[args.chrom_col]
        peak_summit_pos=row[args.peak_start_col]+row[args.summit_col]
        key=(chrom, peak_summit_pos)
        cur_deepshap=foreground_shap_scaled[key]
        tf_start=args.summit_offset+row[args.tf_start_col]-peak_summit_pos
        if tf_start < 0:
            continue 
        tf_end=args.summit_offset+row[args.tf_end_col]-peak_summit_pos
        if tf_end > (2*args.summit_offset):
            continue 
        tf_deepshap=cur_deepshap[tf_start:tf_end]
        tf_name=row[args.tf_name_col]
        vierstra_score=row[args.vierstra_score_col]
        outf.write('\t'.join([str(token) for token in [chrom,
                                                  tf_start,
                                                  tf_end,
                                                  tf_name,
                                                  vierstra_score,
                                                  np.mean(tf_deepshap),
                                                  np.median(tf_deepshap),
                                                  np.percentile(tf_deepshap,90),
                                                  np.max(tf_deepshap),
                                                  np.min(tf_deepshap)]])+'\n')
    outf.close()

if __name__=="__main__":
    main()
    
