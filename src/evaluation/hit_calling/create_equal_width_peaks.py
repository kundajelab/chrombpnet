import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--bed")
    parser.add_argument("--input_len", type=int, default=1000)
    parser.add_argument("--outf")
    return parser.parse_args()

check_dups=[]

def main():

	args=parse_args()
	data=open(args.bed,'r').read().strip().split('\n')
	input_len=args.input_len
	outf=open(args.outf,'w') 

	for region in data:

		tokens=region.split('\t')
		chrom=tokens[0]
		start_pos=int(tokens[1])
		end_pos=int(tokens[2])
		observed_region_size=end_pos-start_pos

		if observed_region_size<=input_len+500:
	        	outf.write(chrom+'\t'+str(start_pos)+'\t'+str(end_pos)+'\n')
		else:
			binned = np.arange(start_pos, end_pos, input_len, dtype=int)
			for idx in range(len(binned)-1):
				if idx < (len(binned)-2):
					outf.write(chrom+'\t'+str(binned[idx])+'\t'+str(binned[idx+1])+'\n')
				else:
					if end_pos - binned[len(binned)-1]<=500:
						outf.write(chrom+'\t'+str(binned[len(binned)-2])+'\t'+str(end_pos)+'\n')
					else:	
						offset1=int((end_pos - binned[len(binned)-1])/2)
						outf.write(chrom+'\t'+str(binned[len(binned)-2])+'\t'+str(binned[len(binned)-2]+500+offset1)+'\n')
						outf.write(chrom+'\t'+str(binned[len(binned)-2]+500+offset1)+'\t'+str(end_pos)+'\n')


#		if observed_region_size<=input_len:
#	        	outf.write(chrom+'\t'+str(start_pos)+'\t'+str(start_pos+input_len)+'\n')
#		else:
#			binned = np.arange(start_pos, end_pos, input_len, dtype=int)

			#offset = binned[-1]+input_len-end_pos
			#if offset>0: # share the offset so that we never fall of the edges on both sides
		#		binned = np.arange(start_pos-offset//2, end_pos-offset//2, input_len, dtype=int)

#			for idx in range(len(binned)):
#				if idx < (len(binned)-1):		
#	        			outf.write(chrom+'\t'+str(binned[idx])+'\t'+str(binned[idx+1])+'\n')
#				else:
#					if end_pos-binned[idx]
#	        			outf.write(chrom+'\t'+str(binned[idx])+'\t'+str(end_pos)+'\n')

	outf.close()
    

if __name__=="__main__":
    main()
