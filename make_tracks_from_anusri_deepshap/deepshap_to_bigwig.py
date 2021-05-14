from os import listdir
from os.path import isfile, join
import operator
import argparse
import pickle
import pyBigWig
import numpy as np

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--path_to_pickle")
    parser.add_argument("--chromsizes",default="/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes")
    parser.add_argument("--output_dir")
    return parser.parse_args()

def main():
    args=parse_args()
    #prepare chromsizes for bigwig header use
    chromsizes=[i.split('\t') for i in open(args.chromsizes,'r').read().strip().split('\n')]
    chromsizes=[tuple([i[0],int(i[1])]) for i in chromsizes]
    chromsizes=sorted(chromsizes,key=operator.itemgetter(0))
    data=pickle.load(open(args.path_to_pickle,'rb'))
    cur_file=args.path_to_pickle.split('/')[-1]
    outf_label_prof=pyBigWig.open('/'.join([args.output_dir,cur_file+'.label_prof.bw']),'w')
    outf_label_prof.addHeader(chromsizes)
    outf_pred_prof=pyBigWig.open('/'.join([args.output_dir,cur_file+'.pred_prof.bw']),'w')
    outf_pred_prof.addHeader(chromsizes)
    outf_profile_shap=pyBigWig.open('/'.join([args.output_dir,cur_file+'.shap_profile.bw']),'w')
    outf_profile_shap.addHeader(chromsizes)
    outf_count_shap=pyBigWig.open('/'.join([args.output_dir,cur_file+'.shap_count.bw']),'w')
    outf_count_shap.addHeader(chromsizes)

    label_prof=data['label_prof']
    pred_prof=data['pred_prof']
    pred_sum=data['pred_sum']
    label_sum=data['label_sum'] 
    seq=data['seq'] 
    profile_shap=data['profile_shap']
    count_shap=data['count_shap']

    print("iterating keys")
    total_keys=len(label_prof)
    key_index=0
    for key in sorted(label_prof):
        print(str(key_index))
        key_index+=1
        chrom=key[0]
        summit=key[1]
        
        label_sum_val=np.exp(np.squeeze(label_sum[key])) # counts
        pred_sum_val=np.exp(np.squeeze(pred_sum[key])) # counts

        label_prof_val=label_prof[key]*label_sum_val
        pred_prof_val=pred_prof[key]*pred_sum_val

        seq_val=seq[key]
        
        if (len(profile_shap[key].shape)>1) and (profile_shap[key].shape[1]==4):
                profile_shap_val=np.sum(profile_shap[key]*seq_val,axis=1)
                count_shap_val=np.sum(count_shap[key]*seq_val,axis=1) 
        else:
            profile_shap_val=profile_shap[key]
            count_shap_val=count_shap[key] 
        n_entries_output=label_prof_val.shape[0]
        start_pos_output=summit-n_entries_output//2
        n_entries_interpretation=profile_shap_val.shape[0]
        start_pos_interpretation=summit-n_entries_interpretation//2
        added=False
        offset=0 
        while not added: 
            try:
                outf_label_prof.addEntries(chrom,start_pos_output+offset,values=label_prof_val[offset::],span=1,step=1)
                outf_pred_prof.addEntries(chrom,start_pos_output+offset,values=pred_prof_val[offset::],span=1,step=1)
                added=True
            except:
                offset+=1
                #print("offset:"+str(offset))
                #print(str(label_prof_val[offset::].shape))
                if offset > label_prof_val.shape[0]:
                    raise Exception()
        added=False
        offset=0
        while not added:
            try:
                outf_profile_shap.addEntries(chrom,start_pos_interpretation+offset,values=profile_shap_val[offset::],span=1,step=1)
                outf_count_shap.addEntries(chrom,start_pos_interpretation+offset,values=count_shap_val[offset::],span=1,step=1)
                added=True
            except:
                offset+=1
                if offset > profile_shap_val.shape[0]:
                    raise Exception()
    outf_label_prof.close()
    outf_pred_prof.close()
    outf_profile_shap.close()
    outf_count_shap.close()
            
                
            

if __name__=="__main__":
    main()
    
