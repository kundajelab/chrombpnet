import argparse
from kerasAC.splits import *
import pandas as pd
import pysam


def parse_args():
    parser=argparse.ArgumentParser(description="generate positive inputs for svm")
    parser.add_argument("--narrowPeak")
    parser.add_argument("--ntrain",type=int,default=None)
    parser.add_argument("--ntest",type=int,default=None,help="None means that all peak in training split will be used")
    parser.add_argument("--out_prefix")
    parser.add_argument("--genome")
    parser.add_argument("--nosort",action='store_true',default=False)
    parser.add_argument("--no_train",action='store_true',default=False)
    parser.add_argument("--no_test", action='store_true',default=False)
    parser.add_argument("--folds",nargs="+",type=int,default=None) 
    return parser.parse_args()

def main():
    args=parse_args()
    #load the narrowPeak file
    peaks=pd.read_csv(args.narrowPeak,header=None,sep='\t')
    #sort by column 9 (8 0-indexed) , which is the qValue, from highest to lowest
    if args.nosort is False:
        sorted_peaks=peaks.sort_values(by=8,ascending=False)
    else:
        sorted_peaks=peaks
    print(sorted_peaks.head())

    #keep dicitionaries to store output file pointers and number of examples processed
    processed=dict()
    outputs=dict()
    fold_dict=dict()
    if args.folds is None: 
        folds=list(splits[args.genome].keys())
    else:
        folds=args.folds
    max_train=args.ntrain
    if max_train is None:
        max_train=sorted_peaks.shape[0]
    max_test=args.ntest
    if max_test is None:
        max_test=sorted_peaks.shape[0] 

    for fold in folds:
        print("setting up for fold:"+str(fold))
        processed[fold]=dict()
        outputs[fold]=dict() 
        if args.no_train is False:
            processed[fold]['train']=0 
            outputs[fold]['train']=open(args.out_prefix+'.train.'+str(fold),'w')
        if args.no_test is False:
            processed[fold]['test']=0
            outputs[fold]['test']=open(args.out_prefix+'.test.'+str(fold),'w')
        args.fold=fold
        args.train_chroms=None
        args.validation_chroms=None
        args.predict_chroms=None
        train_chroms=get_chroms(args,'train')
        valid_chroms=get_chroms(args,'valid')
        test_chroms=get_chroms(args,'test')
        train_chroms=train_chroms+valid_chroms
        fold_dict[fold]={}
        for chrom in train_chroms:
            fold_dict[fold][chrom]='train'
        for chrom in test_chroms:
            fold_dict[fold][chrom]='test'
            
    for index,row in sorted_peaks.iterrows():
        if index%1000==0:
            print(index) 
        cur_chrom=row[0]
        line=cur_chrom+'\t'+str(row[1])+'\t'+str(row[2])+'\n'
        for fold in folds:
            if args.no_train is False:
                if fold_dict[fold][cur_chrom]=='train':
                    if processed[fold]['train']<max_train:
                        processed[fold]['train']+=1
                        #add to train list
                        outputs[fold]['train'].write(line)
            if args.no_test is False:
                if fold_dict[fold][cur_chrom]=='test':
                    if processed[fold]['test']<max_test:
                        processed[fold]['test']+=1
                        outputs[fold]['test'].write(line)
    if args.no_train is False:
        for fold in folds: 
            outputs[fold]['train'].close()
    if args.no_test is False:
        for fold in folds:
            outputs[fold]['test'].close() 
    
if __name__=="__main__":
    main()
    
    
