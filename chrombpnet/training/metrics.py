import numpy as np 
import argparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import jensenshannon
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from chrombpnet.training.utils.metrics_utils import * 

plt.rcParams["figure.figsize"]=10,5
font = {'weight' : 'bold',
        'size'   : 10}
matplotlib.rc('font', **font)

   
def counts_metrics(labels,preds,outf,title):
    '''
    Get count metrics
    '''
    spearman_cor=spearmanr(labels,preds)[0]
    pearson_cor=pearsonr(labels,preds)[0]  
    mse=((labels - preds)**2).mean(axis=0)

    #print("spearman:"+str(spearman_cor))
    #print("pearson:"+str(pearson_cor))
    #print("mse:"+str(mse))

    plt.rcParams["figure.figsize"]=8,8
    fig=plt.figure() 
    density_scatter(labels,
                    preds,
                    xlab='Log Count Labels',
                    ylab='Log Count Predictions')
    plt.suptitle(title+" count: spearman R="+str(round(spearman_cor,3))+", Pearson R="+str(round(pearson_cor,3))+", mse="+str(round(mse,3)))
    plt.legend(loc='best')
    plt.savefig(outf+'.counts_pearsonr.png',format='png',dpi=300)
    
    return spearman_cor, pearson_cor, mse

def profile_metrics(true_counts,pred_probs,pseudocount=0.001):
    '''
    Get profile metrics
    '''
    mnll_pw = []
    mnll_norm = []

    jsd_pw = []
    jsd_norm = []
    jsd_rnd = []
    jsd_rnd_norm = []
    mnll_rnd = []
    mnll_rnd_norm = []

    num_regions = true_counts.shape[0]
    for idx in range(num_regions):
        # mnll
        #curr_mnll = mnll(true_counts[idx,:],  probs=pred_probs[idx,:])
        #mnll_pw.append(curr_mnll)
        # normalized mnll
        #min_mnll, max_mnll = mnll_min_max_bounds(true_counts[idx,:])
        #curr_mnll_norm = get_min_max_normalized_value(curr_mnll, min_mnll, max_mnll)
        #mnll_norm.append(curr_mnll_norm)

        # jsd
        cur_jsd=jensenshannon(true_counts[idx,:]/(pseudocount+np.nansum(true_counts[idx,:])),pred_probs[idx,:])
        jsd_pw.append(cur_jsd)
        # normalized jsd
        min_jsd, max_jsd = jsd_min_max_bounds(true_counts[idx,:])
        curr_jsd_norm = get_min_max_normalized_value(cur_jsd, min_jsd, max_jsd)
        jsd_norm.append(curr_jsd_norm)

        # get random shuffling on labels for a worst case performance on metrics - labels versus shuffled labels
        shuffled_labels=np.random.permutation(true_counts[idx,:])
        shuffled_labels_prob=shuffled_labels/(pseudocount+np.nansum(shuffled_labels))

        # mnll random
        #curr_rnd_mnll = mnll(true_counts[idx,:],  probs=shuffled_labels_prob)
        #mnll_rnd.append(curr_rnd_mnll)
        # normalized mnll random
        #curr_rnd_mnll_norm = get_min_max_normalized_value(curr_rnd_mnll, min_mnll, max_mnll)
        #mnll_rnd_norm.append(curr_rnd_mnll_norm)   

        # jsd random
        curr_jsd_rnd=jensenshannon(true_counts[idx,:]/(pseudocount+np.nansum(true_counts[idx,:])),shuffled_labels_prob)
        jsd_rnd.append(curr_jsd_rnd)
        # normalized jsd random
        curr_rnd_jsd_norm = get_min_max_normalized_value(curr_jsd_rnd, min_jsd, max_jsd)
        jsd_rnd_norm.append(curr_rnd_jsd_norm)

    return np.array(mnll_pw), np.array(mnll_norm), np.array(jsd_pw), np.array(jsd_norm), np.array(jsd_rnd), np.array(jsd_rnd_norm), np.array(mnll_rnd), np.array(mnll_rnd_norm)

def plot_histogram(region_jsd, shuffled_labels_jsd, output_prefix, title):

    #generate histogram distributions 
    num_bins=100
    plt.rcParams["figure.figsize"]=8,8
    
    #plot mnnll histogram 
    #plt.figure()
    #n,bins,patches=plt.hist(mnnll_vals,num_bins,facecolor='blue',alpha=0.5,label="Predicted vs Labels")
    #n1,bins1,patches1=plt.hist(shuffled_labels_mnll,num_bins,facecolor='black',alpha=0.5,label='Shuffled Labels vs Labels')
    #plt.xlabel('Multinomial Negative LL Profile Labels and Predictions in Probability Space')
    #plt.title("MNNLL: "+ tile)
    #plt.legend(loc='best')
    #plt.savefig(output_prefix+".mnnll.png",format='png',dpi=300)
    
    #plot jsd histogram
    plt.figure()
    n,bins,patches=plt.hist(region_jsd,num_bins,facecolor='blue',alpha=0.5,label="Predicted vs Labels")
    n1,bins1,patches1=plt.hist(shuffled_labels_jsd,num_bins,facecolor='black',alpha=0.5,label='Shuffled Labels vs Labels')
    plt.xlabel('Jensen Shannon Distance Profile Labels and Predictions in Probability Space')
    plt.title("JSD Dist: "+title)
    plt.legend(loc='best')
    plt.savefig(output_prefix+".profile_jsd.png",format='png',dpi=300)


