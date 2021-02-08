import pandas as pd
import random

regions=open('regions_to_score.bed','r').read().strip().split('\n')
outf=open("regions.bed",'w')
for line in regions:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=int(tokens[1])
    summit=start+int(tokens[-2])
    start=summit-500
    end=summit+500
    outf.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+chrom+":"+str(summit)+'\n')
outf.close()

    
