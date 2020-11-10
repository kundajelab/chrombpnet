import pickle
import numpy as np 
data=pickle.load(open("K562.DNASE.bias_corrected_bpnet_tobias.unplugged.fold0.deepSHAP",'rb'))
outf=open('K562.DNASE.DeepSHAP.fold0.tsv','w')
outf.write('Chrom\tPos\tProfileDeepSHAP\tCountDeepSHAP\n')
keys=list(data['profile_shap'].keys())
count=0
for key in keys:
    count+=1
    print(str(count))
    profile_score=data['profile_shap'][key]
    count_score=data['count_shap'][key]
    seq=data['seq'][key]
    profile_score_scaled=np.sum(profile_score*seq,axis=1)
    count_score_scaled=np.sum(count_score*seq,axis=1)
    chrom=key[0]
    summit=key[1]
    start_pos=summit-500
    for i in range(1000):
        outf.write(chrom+'\t'+str(start_pos+i)+'\t'+str(profile_score_scaled[i])+'\t'+str(count_score_scaled[i])+'\n')
outf.close()
