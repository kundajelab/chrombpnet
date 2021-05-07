from kerasAC.util import *
import argparse
import Bio
from Bio import motifs
import numpy as np
import pickle
import pdb

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--meme_file",default="JASPAR2018_CORE_vertebrates_non-redundant_pfms.meme")
    parser.add_argument("--top_n",default=3,type=int)
    parser.add_argument("--out_pickle",default='JASPER.top3.pickle')
    return parser.parse_args()

def main():
    args=parse_args()
    pos_to_letter={0:'A',
                   1:'C',
                   2:'G',
                   3:'T'}
    letter_to_pos={'A':0,
                   'C':1,
                   'G':2,
                   'T':3}
    with open(args.meme_file) as f:
        pwms=motifs.parse(f,'MINIMAL')
        
    pwm_dict={}
    for pwm in pwms:
        vals=np.array([list(pwm.pwm['A']),
                       list(pwm.pwm['C']),
                       list(pwm.pwm['G']),
                       list(pwm.pwm['T'])])
        consensus=pwm.consensus
        name=pwm.name 
        seq_to_prob={}
        consensus=list(consensus)
        probs=np.max(vals,axis=0)
        for j in range(vals.shape[1]):
            for i in range(4):
                new_letter=pos_to_letter[i]
                new_consensus=consensus
                new_consensus[j]=new_letter
                new_probs=probs 
                new_probs[j]==vals[i,j]
                prob_observed=np.sum([np.log(z) for z in new_probs])
                seq_to_prob[''.join(new_consensus)]=prob_observed 
        #sort seq_to_prob by prob, in descending order 
        most_likely_seqs=sorted(seq_to_prob.items(), key=lambda item: item[1],reverse=True)[0:args.top_n]
        print(most_likely_seqs)
        pwm_dict[name]=[seq[0] for seq in most_likely_seqs]
    #pickle the pwm dict 
    with open(args.out_pickle,'wb') as f:
        pickle.dump(pwm_dict,f,pickle.HIGHEST_PROTOCOL)
        
        

if __name__=='__main__':
    main()
    
