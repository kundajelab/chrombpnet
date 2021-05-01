#!/usr/bin/env python
# coding: utf-8

# In[39]:
import matplotlib 
from matplotlib import pyplot as plt
import pandas as pd 
import pickle
import pybedtools 
import numpy as np 
from scipy.stats import gmean


# In[40]:

import sys 
cell=sys.argv[1]


# In[2]:


background_hits=pickle.load(open("../interpret/"+cell+".DNASE.bias_corrected_bpnet_tobias.unplugged.SHUFFLED.fold0.deepSHAP",'rb'))


# In[3]:


foreground_hits=pickle.load(open("../interpret/"+cell+".DNASE.bias_corrected_bpnet_tobias.unplugged.fold0.deepSHAP",'rb'))


# In[4]:


foreground_shap=foreground_hits['profile_shap']
foreground_seq=foreground_hits['seq']
foreground_shap_scaled={} 
for key in foreground_shap: 
    foreground_shap_scaled[key]=np.sum(foreground_shap[key]*foreground_seq[key],axis=1)


# In[5]:


background_shap=background_hits['profile_shap']
keys=list(background_shap.keys())
import random 
random.seed(1234)
key_subset=random.sample(keys,10000)
background_flat=[] 
for key in key_subset: 
    background_flat=background_flat+background_shap[key].tolist()


# In[9]:


background_flat=np.asarray(background_flat)


# In[11]:


#restrict to only positive scores 
background_flat=background_flat[background_flat>=0]


# In[15]:


#plot the background distribution 
plt.hist(background_flat,bins=1000)
plt.yscale("log")
plt.title(cell+" shuffled background profile deepSHAP")
plt.savefig(cell+".shuffled.background.profile.deepshap.png")


# In[16]:


from statsmodels.distributions.empirical_distribution import ECDF


# In[17]:


ecdf=ECDF(abs(background_flat))


# In[18]:

plt.figure()
plt.plot(ecdf.x, ecdf.y)
plt.title(cell+" shuffled background profile deepSHAP CDF")
plt.savefig(cell+".shuffled.background.profile.deepshap.cdf.png")


# In[20]:


#pickle the background scores 
import pickle 
# open a file, where you ant to store the data
file = open(cell+'.dnase.background.deepshap.ecdf.pickle', 'wb')
# dump information to that file
pickle.dump(ecdf, file)
file.close()


# In[36]:


intersections=pd.read_csv(cell+'.idr.peaks.50bp.around.summit.overlap.tf.all.bed',header=None,sep='\t')


# In[37]:


intersections.head()


# In[ ]:


import pdb
tf_to_geom_mean_pval={} 
tf_to_90_percentile_pval={} 
tf_to_pval_of_mean={} 
outf=open(cell+".tf.sig.tsv",'w')
count=0 
for key in foreground_shap_scaled: 
    if count %100==0: 
        print(str(count))
    count+=1
    chrom=key[0]
    summit=key[1] 
    deepshap_vals=foreground_shap_scaled[key]
    #change negative values to 0; not associated w/ predictive instances  
    deepshap_vals[deepshap_vals<0]=0 
    #get empirical p-values 
    deepshap_pvals=1-ecdf(abs(deepshap_vals))

    hits=intersections[(intersections[7]==chrom) & (intersections[10]==summit)]
    num_hits=hits.shape[0]
    tf_names=hits[3]
    tf_start=1056+hits[1]-hits[10]
    tf_end=1056+hits[2]-hits[10]
    for i in range(num_hits): 
        cur_tf=tf_names.iloc[i]
        cur_tf_start=tf_start.iloc[i] 
        cur_tf_end=tf_end.iloc[i]
        if cur_tf not in tf_to_geom_mean_pval: 
            tf_to_geom_mean_pval[cur_tf]=[]
            tf_to_90_percentile_pval[cur_tf]=[]
            tf_to_pval_of_mean[cur_tf]=[] 
            
        #pvalue of the mean deepSHAP score,make sure pval >=1e-320 to avoid numerical precision issues with log 
        pval_of_mean=max([1-ecdf(np.mean(deepshap_vals[cur_tf_start:cur_tf_end])),1e-320])
        geom_mean_pval=max([gmean(deepshap_pvals[cur_tf_start:cur_tf_end]),1e-320])
        pval_90_percentile=max([np.percentile(deepshap_pvals[cur_tf_start:cur_tf_end],90),1e-320])
        
        #compute -log10(pval)
        pval_of_mean=-1*np.log10(pval_of_mean)
        geom_mean_pval=-1*np.log10(geom_mean_pval)
        pval_90_percentile=-1*np.log10(pval_90_percentile)
        
        
        #append to list for distribution 
        tf_to_pval_of_mean[cur_tf].append(pval_of_mean)
        tf_to_geom_mean_pval[cur_tf].append(geom_mean_pval)
        tf_to_90_percentile_pval[cur_tf].append(pval_90_percentile)
        #keep track of significant hits
        if max([pval_of_mean,geom_mean_pval,pval_90_percentile]) > 2: 
            outf.write(chrom+'\t'+str(hits[1].iloc[i])+'\t'+str(hits[2].iloc[i])+'\t'+str(hits[5].iloc[i])+'\t'+cur_tf+'\t'+str(round(pval_of_mean,2))+'\t'+str(round(geom_mean_pval,2))+'\t'+str(round(pval_90_percentile,2))+'\n')
outf.close()


# In[ ]:


import pickle
with open(cell+'.tfscan.tf_to_geom_mean_pval.pickle', 'wb') as handle:
    pickle.dump(tf_to_geom_mean_pval, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(cell+'.tfscan.tf_to_90_percentile_pval.pickle', 'wb') as handle:
    pickle.dump(tf_to_90_percentile_pval, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(cell+'.tfscan.tf_to_pval_of_mean.pickle', 'wb') as handle:
    pickle.dump(tf_to_pval_of_mean, handle, protocol=pickle.HIGHEST_PROTOCOL)

keys=list(tf_to_geom_mean_pval.keys())
keys.sort()

plt.rcParams["figure.figsize"]=5,10
for key in keys:
    print(str(key))
    fig, axes = plt.subplots(3, 1)
    cur_to_geom_mean_pval=tf_to_geom_mean_pval[key]
    cur_to_90_percentile_pval=tf_to_90_percentile_pval[key]
    cur_to_pval_of_mean=tf_to_pval_of_mean[key]
    axes[0].hist(cur_to_geom_mean_pval,bins=150)
    axes[0].set_title("Geometric mean -log10(Pval)")
    axes[1].hist(cur_to_90_percentile_pval,bins=150)
    axes[1].set_title("90% DeepSHAP -log10(Pval)")
    axes[2].hist(cur_to_pval_of_mean,bins=150) 
    axes[2].set_title("-log10(Pval of mean deepSHAP)")
    axes[0].set_xlim(0,10)
    axes[1].set_xlim(0,10)
    axes[2].set_xlim(0,10)
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[2].set_yscale('log') 
    plt.suptitle(cell+' DNASE, '+str(key))
    plt.savefig(key.replace('/','-')+'.'+cell+'.idr.png')
