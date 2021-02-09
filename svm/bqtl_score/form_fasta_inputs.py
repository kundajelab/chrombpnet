import argparse
import pysam
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--snp_file")
    parser.add_argument("--chrom_col")
    parser.add_argument("--pos_col")
    parser.add_argument("--ref_col")
    parser.add_argument("--alt_col")
    parser.add_argument("--ref_fasta")
    parser.add_argument("--flank_size",type=int,default=500)
    parser.add_argument("--offset_1",action="store_true",default=False)
    parser.add_argument("--rsid_col",help="column that stores rsid or similar")
    parser.add_argument("--out_prefix") 
    return parser.parse_args()

def main():
    args=parse_args()
    ref=pysam.FastaFile(args.ref_fasta)
    snps=pd.read_csv(args.snp_file,sep='\t',header=0)
    out_ref=open(args.out_prefix+".ref.fa",'w')
    out_alt=open(args.out_prefix+".alt.fa",'w')
    for index,row in snps.iterrows():
        cur_chrom=row[args.chrom_col]
        cur_pos=row[args.pos_col]
        ref_allele=row[args.ref_col]
        alt_allele=row[args.alt_col]
        snp_id=row[args.rsid_col]
        if args.offset_1 is True:
            #shift position by 1
            cur_pos-=1
        total_length=2*args.flank_size
        left_flank=ref.fetch(cur_chrom,cur_pos-args.flank_size,cur_pos)
        right_flank=ref.fetch(cur_chrom,cur_pos+1,cur_pos+args.flank_size)
        if len(ref_allele)>1:
            print("warning:the ref allele for rsid:"+args.snp_id+":"+ref_allele)
            ref_allele=ref_allele[0]
        if len(alt_allele)>1:
            print("warning:the alt allele for rsid:"+args.snp_id+":"+alt_allele)
            alt_allele=alt_allele[0]
        ref_seq=left_flank+ref_allele+right_flank
        alt_seq=left_flank+alt_allele+right_flank
        out_ref.write(">"+snp_id+":ref:"+ref_allele+'\n'+ref_seq+'\n')
        out_alt.write(">"+snp_id+":alt:"+alt_allele+'\n'+alt_seq+'\n')
    out_ref.close()
    out_alt.close()

if __name__=="__main__":
    main()
    
    
