import pdb 
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="aggregate jsd/mnnll summaries")
    parser.add_argument("--prefix",nargs="+")
    parser.add_argument("--files",nargs="+")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args() 
    jsd_dict=dict()
    mnnll_dict=dict()
    score_sets=set([])
    for i in range(len(args.files)):
        cur_prefix=args.prefix[i]
        cur_file=open(args.files[i],'r').read().strip().split('\n')[1::]
        print(str(args.files[i]))
        score_sets.add(cur_prefix)
        for line in cur_file:
            tokens=line.split('\t')
            if len(tokens)!=5:
                print(line)
                continue
            pos=tuple(tokens[0:2])
            jsd=tokens[2]
            mnnll=tokens[4]
            if pos not in jsd_dict:
                jsd_dict[pos]={}
            if pos not in mnnll_dict:
                mnnll_dict[pos]={}
            jsd_dict[pos][cur_prefix]=jsd
            mnnll_dict[pos][cur_prefix]=mnnll
    score_sets=list(score_sets)

    outf_jsd=open(args.outf+".jsd.txt",'w')
    outf_jsd.write('chrom\tsummit\t'+'\t'.join(score_sets)+'\n')
    for pos in jsd_dict:
        if len(jsd_dict[pos].keys())!=len(score_sets):
            continue
        outf_jsd.write(pos[0]+'\t'+pos[1])
        for score_set in score_sets:
            outf_jsd.write('\t'+jsd_dict[pos][score_set])
        outf_jsd.write('\n')
    outf_jsd.close()
    
    outf_mnll=open(args.outf+".mnnll.txt",'w')
    outf_mnll.write('chrom\tsummit\t'+'\t'.join(score_sets)+'\n')
    for pos in mnnll_dict:
        if len(mnnll_dict[pos].keys())!=len(score_sets):
            continue
        outf_mnll.write(pos[0]+'\t'+pos[1])
        for score_set in score_sets:
            outf_mnll.write('\t'+mnnll_dict[pos][score_set])
        outf_mnll.write('\n')
    outf_mnll.close()
    

if __name__=="__main__":
    main()
    
