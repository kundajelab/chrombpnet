import sys
fold=sys.argv[1]
data=open('labels.'+fold+'.bed','r').read().strip().split('\n')
outf=open('labels.'+fold+'.filled.bed','w')
for line in data:
    tokens=line.split('\t')
    if tokens[-1]!="":
        outf.write(line+'\n')
    else:
        outf.write(line+'NA\n')
outf.close()
