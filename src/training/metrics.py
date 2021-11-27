import numpy as np 
import argparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import jensenshannon
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from utils.metrics_utils import * 

plt.rcParams["figure.figsize"]=10,5
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}
matplotlib.rc('font', **font)

   
def counts_metrics(labels,preds,outf,title):
    '''
    Get count metrics
    '''
    spearman_cor=spearmanr(labels,preds)[0]
    pearson_cor=pearsonr(labels,preds)[0]  
    mse=((labels - preds)**2).mean(axis=0)

    print("spearman:"+str(spearman_cor))
    print("pearson:"+str(pearson_cor))
    print("mse:"+str(mse))

    plt.rcParams["figure.figsize"]=8,8
    fig=plt.figure() 
    density_scatter(labels,
                    preds,
                    xlab='Log Count Labels',
                    ylab='Log Count Predictions')
    plt.suptitle(title+" count: spearman R="+str(round(spearman_cor,3))+", Pearson R="+str(round(pearson_cor,3))+", mse="+str(round(mse,3)))
    plt.legend(loc='best')
    plt.savefig(outf+'.png',format='png',dpi=300)
    
    return spearman_cor, pearson_cor, mse

def profile_metrics(true_counts,pred_probs):
    '''
    Get profile metrics
    '''
    mnll_pw = []
    mnll_norm = []

    jsd_pw = []
    jsd_norm = []
    jsd_rnd = []
    jsd_rnd_norm = []

    num_regions = true_counts.shape[0]
    for idx in range(num_regions):
        curr_mnll = mnll(true_counts[idx,:],  probs=pred_probs[idx,:])
        min_mnll, max_mnll = mnll_min_max_bounds(true_counts[idx,:])
        curr_mnll_norm = get_min_max_normalized_value(curr_mnll, min_mnll, max_mnll)
        mnll_pw.append(curr_mnll)
        mnll_norm.append(curr_mnll_norm)

        cur_jsd=jensenshannon(true_counts[idx,:]/np.nansum(true_counts[idx,:]),pred_probs[idx,:])
        min_jsd, max_jsd = jsd_min_max_bounds(true_counts[idx,:])
        curr_jsd_norm = get_min_max_normalized_value(cur_jsd, min_jsd, max_jsd)
        jsd_pw.append(cur_jsd)
        jsd_norm.append(curr_jsd_norm)

        shuffled_labels=np.random.permutation(true_counts[idx,:])
        shuffled_labels_prob=shuffled_labels/np.nansum(shuffled_labels)
        curr_jsd_rnd=jensenshannon(true_counts[idx,:]/np.nansum(true_counts[idx,:]),shuffled_labels_prob)
        jsd_rnd.append(curr_jsd_rnd)
        curr_rnd_jsd_norm = get_min_max_normalized_value(curr_jsd_rnd, min_jsd, max_jsd)
        jsd_rnd_norm.append(curr_rnd_jsd_norm)

    return np.array(mnll_pw), np.array(mnll_norm), np.array(jsd_pw), np.array(jsd_norm), np.array(jsd_rnd), np.array(jsd_rnd_norm)


