import pickle
import numpy as np
import pyBigWig
import operator
data=pickle.load(open("hg19.K562.DNASE.bias_corrected_bpnet_tobias.unplugged.fold0.deepSHAP",'rb'))
print("loaded pickle of deepshap scores")

bw_profile=pyBigWig.open('K562.DNASE.DeepSHAP.fold0.hg19.idr.profile.bw','w')
bw_count=pyBigWig.open('K562.DNASE.DeepSHAP.fold0.hg19.idr.count.bw','w')
chromsizes=[i.split('\t') for i in open("hg19.chrom.sizes",'r').read().strip().split('\n')]
chromsizes=[tuple([i[0],int(i[1])]) for i in chromsizes]
chromsizes=sorted(chromsizes,key=operator.itemgetter(0))
bw_profile.addHeader(chromsizes)
bw_count.addHeader(chromsizes)

keys=list(data['profile_shap'].keys())
keys.sort()
count=0
last_end=0
last_chrom=0 
for key in keys:
    count+=1
    print(str(count))
    profile_score=data['profile_shap'][key]
    count_score=data['count_shap'][key]
    seq=data['seq'][key]
    profile_score_scaled=np.sum(profile_score*seq,axis=1)
    count_score_scaled=np.sum(count_score*seq,axis=1)
    cur_chrom=key[0]
    summit=key[1]
    cur_start=summit-1057
    if (last_chrom==cur_chrom) and (cur_start < last_end):
        delta=last_end-cur_start
        profile_score_scaled=profile_score_scaled[delta::]
        count_score_scaled=count_score_scaled[delta::]        
        cur_start=last_end 
    cur_end=summit+1057
    last_end=cur_end
    last_chrom=cur_chrom
    if len(profile_score_scaled)==0:
        continue
    print("cur_start:"+str(cur_start))
    print("cur_end:"+str(cur_end))
    print(str(len(profile_score_scaled)))
    
    bw_profile.addEntries(cur_chrom,cur_start,values=profile_score_scaled,span=1,step=1)
    bw_count.addEntries(cur_chrom,cur_start,values=count_score_scaled,span=1,step=1)
bw_profile.close()
bw_count.close()
