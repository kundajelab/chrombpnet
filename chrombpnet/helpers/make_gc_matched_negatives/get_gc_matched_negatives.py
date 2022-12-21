import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import random
import csv 
import json
import sys

def parse_args():
    parser=argparse.ArgumentParser(description="generate a bed file of non-peak regions that are gc-matched with foreground")
    parser.add_argument("-c","--candidate_negatives",help="candidate negatives bed file with gc content in 4th column rounded to 2 decimals")
    parser.add_argument("-f","--foreground_gc_bed", help="regions with their corresponding gc fractions for matching, 4th column has gc content value rounded to 2 decimals")
    parser.add_argument("-o","--output_prefix", help="gc-matched non-peaks output file name")
    parser.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
    parser.add_argument("-npr", "--neg_to_pos_ratio_train", type=int, default=1, help="Ratio of negatives to positives to sample for training")
    parser.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for generating nonpeaks")

    return parser.parse_args()

def remap_chrom(chrom, splits_dict):
    '''
    Remapping chromosome names - we will not differentiate between the train/valid/tes chromsomes
    when sampling negatives.
    '''
    if chrom in splits_dict["train"]:
        chrom_mod = "chrom_train"
    elif chrom in splits_dict["valid"]:
        chrom_mod = "chrom_valid"
    elif chrom in splits_dict["test"]:
        chrom_mod = "chrom_test"
    else:
        chrom_mod = "ignore"
    return chrom_mod


def make_gc_dict(candidate_negatives, splits_dict):
    """
    Imports the candidate negatives into a dictionary structure.
    The `key` is the gc content fraction, and the `values` are a list 
    containing the (chrom,start,end) of a region with the corresponding 
    gc content fraction.
    """
    data=open(candidate_negatives,'r').readlines()
    gc_dict={}
    index=0 
    ignored_chroms = []
    for line in tqdm(list(data)):
        line=line.strip('\n') 
        index+=1
        tokens=line.split('\t')
        chrom=tokens[0]
        gc=float(tokens[-1])
        start=tokens[1]
        end=tokens[2]
        chrom_real=chrom
        chrom = remap_chrom(chrom, splits_dict)
        if chrom == "ignore":
            ignored_chroms.append(chrom_real)
            continue
        if chrom not in gc_dict:
            gc_dict[chrom]={}
        if gc not in gc_dict[chrom]:
            gc_dict[chrom][gc]=[(chrom,start,end,chrom_real)]
        else:
            gc_dict[chrom][gc].append((chrom,start,end,chrom_real))

    print("Following background chromosomes {}  were ignored since they are not present in the given fold".format(",".join(list(set(ignored_chroms)))))
    return gc_dict

def scale_gc(cur_gc):
    """
    Randomly increase/decrease the gc-fraction value by 0.01
    """
    if random.random()>0.5:
        cur_gc+=0.01
    else:
        cur_gc-=0.01
    cur_gc=round(cur_gc,2)
    if cur_gc<=0:
        cur_gc+=0.01
    if cur_gc>=1:
        cur_gc-=0.01
    assert cur_gc >=0
    assert cur_gc <=1
    return cur_gc 

def adjust_gc(chrom,cur_gc,negatives,used_negatives):
    """
    Function that checks if (1) the given gc fraction value is available
    in the negative candidates or (2) if the given gc fraction value has 
    candidates not already sampled. If eitheir of the condition fails we  
    sample the neighbouring gc_fraction value by randomly scaling with 0.01.
    """
    if chrom  not in used_negatives:
        used_negatives[chrom]={}

    if cur_gc not in used_negatives[chrom]:
        used_negatives[chrom][cur_gc]=[]

    while (cur_gc not in negatives[chrom]) or (len(used_negatives[chrom][cur_gc])>=len(negatives[chrom][cur_gc])):
        cur_gc=scale_gc(cur_gc)
        if cur_gc not in used_negatives[chrom]:
            used_negatives[chrom][cur_gc]=[]
    return cur_gc,used_negatives 

        
    
def main(args): 


    splits_dict=json.load(open(args.chr_fold_path))
    
    random.seed(args.seed)

    negatives=make_gc_dict(args.candidate_negatives, splits_dict)
    used_negatives=dict()
    cur_peaks=pd.read_csv(args.foreground_gc_bed,header=None,sep='\t')
    negatives_bed = []
    print(len(list(cur_peaks.iterrows())))
    
    foreground_gc_vals = []
    output_gc_vals = []
    ignored_chroms = []
    for index,row in tqdm(list(cur_peaks.iterrows())): 

        chrom=row[0]
        start=row[1]
        end=row[2]
        gc_value=row[3]

        chrom_real=chrom
        chrom = remap_chrom(chrom, splits_dict)
        if chrom == "ignore":
            ignored_chroms.append(chrom_real)
            continue

        if chrom=="chrom_train" or chrom=="chrom_valid":
            neg_to_pos_ratio = args.neg_to_pos_ratio_train
        else:
            neg_to_pos_ratio = 1

        # for every gc value in positive how many negatives to find
        # we will keep the ratio of positives to negatives in the test set same
        for rep in range(neg_to_pos_ratio):
            cur_gc,used_negatives=adjust_gc(chrom,gc_value,negatives,used_negatives)
            num_candidates=len(negatives[chrom][cur_gc])
            rand_neg_index=random.randint(0,num_candidates-1)
            while rand_neg_index in used_negatives[chrom][cur_gc]:
                cur_gc,used_negatives=adjust_gc(chrom,cur_gc,negatives,used_negatives)
                num_candidates=len(negatives[chrom][cur_gc])
                rand_neg_index=random.randint(0,num_candidates-1)

            used_negatives[chrom][cur_gc].append(rand_neg_index)
            neg_tuple=negatives[chrom][cur_gc][rand_neg_index]
            neg_chrom=neg_tuple[0]
            neg_start=neg_tuple[1]
            neg_end=neg_tuple[2]
            neg_chrom_real=neg_tuple[3]
            negatives_bed.append([neg_chrom_real,int(neg_start),int(neg_end), cur_gc]) 
            output_gc_vals.append(cur_gc)
            foreground_gc_vals.append(gc_value)       
  
    print("Following foreground chromosomes {} were ignored since they are not present in the given fold".format(",".join(list(set(ignored_chroms)))))     
    negatives_bed = pd.DataFrame(negatives_bed)
    negatives_bed.to_csv(args.output_prefix+".bed", sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    # checking how far the true distribution of foreground is compared to the backgrounds generated
    bins = np.linspace(0, 1, 100)
    plt.hist([output_gc_vals,foreground_gc_vals], bins, density=True, label=['negatives gc distribution', "foreground gc distribution"])
    plt.xlabel("GC content")
    plt.ylabel("Density")
    plt.legend(loc='upper right')
    plt.savefig(args.output_prefix+"_compared_with_foreground.png")

    
if __name__=="__main__":
    args=parse_args()
    main(args)
