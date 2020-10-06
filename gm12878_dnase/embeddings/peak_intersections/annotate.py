peaks_to_tf={}
for f in ['CTCF.bed','STAT1.bed','JUND.bed','NFKB.bed','POU2F2.bed','SPI1.bed']:
    tf=f.split('.')[0]
    data=open(f,'r').read().strip().split('\n')
    for line in data:
        if line not in peaks_to_tf:
            peaks_to_tf[line]=tf
        else:
            peaks_to_tf[line]=','.join([peaks_to_tf[line],tf])
outf=open('peaks_to_tf.bed','w')
for line in peaks_to_tf:
    outf.write(line+'\t'+peaks_to_tf[line]+'\n')
outf.close()
