from os import listdir
from os.path import isfile, join
import operator
import argparse
import pickle
import pyBigWig
import numpy as np
import pandas as pd
from scipy.special import softmax,logit
import os


def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--path_to_pickle")
    parser.add_argument("--chromsizes",default="/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes")
    parser.add_argument("--output_dir")
    return parser.parse_args()

def main():
    args=parse_args()
    #prepare chromsizes for bigwig header use
    chromsizes=[i.split('\t') for i in open(args.chromsizes,'r').read().strip().split('\n')]
    chromsizes=[tuple([i[0],int(i[1])]) for i in chromsizes]
    chromsizes=sorted(chromsizes,key=operator.itemgetter(0))

    profile_preds=pickle.load(open(os.path.join(args.path_to_pickle, "profile_predictions.pkl"),'rb'))
    count_shap=pickle.load(open(os.path.join(args.path_to_pickle, "count_shap.pkl"),'rb'))
    profile_shap=pickle.load(open(os.path.join(args.path_to_pickle, "profile_shap.pkl"),'rb'))
    count_preds = pd.read_csv(os.path.join(args.path_to_pickle,"count_predictions_alt_and_ref.tsv"), header=0,index_col=0,sep='\t')

    ref_profile_bw=pyBigWig.open(os.path.join(args.output_dir,'profile_ref.bw'),'w')
    ref_profile_bw.addHeader(chromsizes)
    alt_profile_bw=pyBigWig.open(os.path.join(args.output_dir,'profile_alt.bw'),'w')
    alt_profile_bw.addHeader(chromsizes)
    diff_profile_bw=pyBigWig.open(os.path.join(args.output_dir,'profile_diff.bw'),'w')
    diff_profile_bw.addHeader(chromsizes)
    ref_profile_shap=pyBigWig.open(os.path.join(args.output_dir,'profile_shap_ref.bw'),'w')
    ref_profile_shap.addHeader(chromsizes)
    alt_profile_shap=pyBigWig.open(os.path.join(args.output_dir,'profile_shap_alt.bw'),'w')
    alt_profile_shap.addHeader(chromsizes)
    ref_count_shap=pyBigWig.open(os.path.join(args.output_dir,'count_shap_ref.bw'),'w')
    ref_count_shap.addHeader(chromsizes)
    alt_count_shap=pyBigWig.open(os.path.join(args.output_dir,'count_shap_alt.bw'),'w')
    alt_count_shap.addHeader(chromsizes)
    
    prev_chr=""
    keys = count_preds["rsid"].values
    myList = sorted([(int(key.split("_")[0].replace("chr","").replace("X","23").replace("Y", "24")), int(key.split("_")[1])) for key in keys])
    sorted_idx = sorted(range(len(myList)),key=myList.__getitem__)

    for key_idx in sorted_idx:
        key = keys[key_idx]
        print(key)
        temp=key.split("_")
        chrom=temp[0]
        Pos0=int(temp[1])
        prof_flank=500
        shap_flank=500
        curr_chrom=chrom
        if prev_chr != curr_chrom:
            cur_pos=0
            prev_chr=curr_chrom

        ref_profile_prof = np.array(softmax(profile_preds[key]['ref'][:,0],axis=0))
        ref_counts = count_preds.loc[key, 'ref']
        ref_profile = ref_profile_prof*np.expand_dims(np.exp([ref_counts]), axis=1)[0]

        alt_profile_prof = np.array(softmax(profile_preds[key]['alt'][:,0],axis=0))
        alt_counts = count_preds.loc[key, 'alt']
        alt_profile = alt_profile_prof*np.expand_dims(np.exp([alt_counts]), axis=1)[0]

        diff_prof_ref_alt = np.log2(ref_profile+1e-8)-np.log2(alt_profile+1e-8)

        profile_shap_ref = profile_shap[key]["ref"]
        profile_shap_alt = profile_shap[key]["alt"]
        count_shap_ref = count_shap[key]["ref"]
        count_shap_alt = count_shap[key]["alt"]

        #print(ref_profile.shape)
        #print(diff_prof_ref_alt.shape)
        #print(profile_shap_ref.shape)
        #print(cur_pos)
        #print(Pos0-int(prof_flank))
        if Pos0-int(prof_flank) >= cur_pos:
            ref_profile_bw.addEntries(chrom,Pos0-int(prof_flank),values=ref_profile,span=1,step=1)
            alt_profile_bw.addEntries(chrom,Pos0-int(prof_flank),values=alt_profile,span=1,step=1)
            diff_profile_bw.addEntries(chrom,Pos0-int(prof_flank),values=diff_prof_ref_alt,span=1,step=1)
            ref_profile_shap.addEntries(chrom,Pos0-int(shap_flank),values=profile_shap_ref[1057-shap_flank:1057+shap_flank],span=1,step=1)
            alt_profile_shap.addEntries(chrom,Pos0-int(shap_flank),values=profile_shap_alt[1057-shap_flank:1057+shap_flank],span=1,step=1)
            ref_count_shap.addEntries(chrom,Pos0-int(shap_flank),values=count_shap_ref[1057-shap_flank:1057+shap_flank],span=1,step=1)
            alt_count_shap.addEntries(chrom,Pos0-int(shap_flank),values=count_shap_alt[1057-shap_flank:1057+shap_flank],span=1,step=1)
            cur_pos = Pos0 + int(prof_flank)
        else:
            offset=cur_pos-(Pos0-int(prof_flank))
            print(offset)
            ref_profile_bw.addEntries(chrom,Pos0-int(prof_flank)+offset,values=ref_profile[offset:],span=1,step=1)
            alt_profile_bw.addEntries(chrom,Pos0-int(prof_flank)+offset,values=alt_profile[offset:],span=1,step=1)
            diff_profile_bw.addEntries(chrom,Pos0-int(prof_flank)+offset,values=diff_prof_ref_alt[offset:],span=1,step=1)
            ref_profile_shap.addEntries(chrom,Pos0-int(shap_flank)+offset,values=profile_shap_ref[1057-shap_flank+offset:1057+shap_flank],span=1,step=1)
            alt_profile_shap.addEntries(chrom,Pos0-int(shap_flank)+offset,values=profile_shap_alt[1057-shap_flank+offset:1057+shap_flank],span=1,step=1)
            ref_count_shap.addEntries(chrom,Pos0-int(shap_flank)+offset,values=count_shap_ref[1057-shap_flank+offset:1057+shap_flank],span=1,step=1)
            alt_count_shap.addEntries(chrom,Pos0-int(shap_flank)+offset,values=count_shap_alt[1057-shap_flank+offset:1057+shap_flank],span=1,step=1)
            cur_pos = Pos0 + int(prof_flank)

    ref_profile_bw.close()
    alt_profile_bw.close()
    diff_profile_bw.close()
    ref_profile_shap.close()
    alt_profile_shap.close()
    ref_count_shap.close()
    alt_count_shap.close()










                
            

if __name__=="__main__":
    main()
