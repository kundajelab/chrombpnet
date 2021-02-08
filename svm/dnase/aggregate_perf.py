import argparse
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--perf_dirs",nargs="+")
    parser.add_argument("--outf",default="aggregate_performance.tsv")
    return parser.parse_args()

def main():
    args=parse_args()
    perf_dict=dict()
    outf=open(args.outf,'w')
    outf.write('Task\tFold\tModel\tTrainSet\tTestSet\tauPRC\n')
    for perf_dir in args.perf_dirs:
        perf_name_tokens=perf_dir.split('.')
        cur_model=perf_name_tokens[1]
        cur_trainset=perf_name_tokens[2]
        cur_testset=perf_name_tokens[3]
        print("model:"+cur_model)
        print("trainset:"+cur_trainset)
        print("testset:"+cur_testset) 
        perf_data=open(perf_dir+'/perf.metrics.txt','r').read().strip().split('\n')
        for line in perf_data[1::]:
            tokens=line.split('\t')
            print(tokens)
            cur_task=tokens[0]
            cur_fold=tokens[1]                                
            cur_auprc=tokens[2]
            outf.write(cur_task+'\t'+cur_fold+'\t'+cur_model+'\t'+cur_trainset+'\t'+cur_testset+'\t'+cur_auprc+'\n')
            
    outf.close()
    
    
if __name__=="__main__":
    main()
    
    
