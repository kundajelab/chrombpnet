import sys
input=open(sys.argv[1],'r').read().strip().split('\n')
#seen=dict()
outf=open(sys.argv[1]+'.bed','w')
for line in input:
    tokens=line.split('\t')
    pos=tokens[0]
    #if pos in seen:
    #    continue
    #seen[pos]=1
    splits1=pos.split(':')
    chrom=splits1[0]
    splits2=splits1[1].split('-')
    val=tokens[1]
    outf.write(chrom+'\t'+splits2[0]+'\t'+splits2[1]+'\t'+val+'\n')
outf.close()

    
