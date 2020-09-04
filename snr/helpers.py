#helper functions for comparing Vierstra footprints to deepSHAP and gkmexplain scores
## for plotting 
import matplotlib 
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"]=10,5
plt.rcParams['axes.xmargin'] = 0

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}

matplotlib.rc('font', **font)
import numpy as np
import pandas as pd
from kerasAC.vis import *
from sklearn.metrics import average_precision_score
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn


def density_scatter(x, y, xlab, ylab, ax = None, sort = True, bins = 20):
    """
    Scatter plot colored by 2d histogram
    """
    bad_indices=np.where(np.isnan(x))+np.where(np.isnan(y))
    x=x[~np.isin(np.arange(x.size),bad_indices)]
    y=y[~np.isin(np.arange(y.size),bad_indices)]

    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0
    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z,s=0.3 )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    #cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    #cbar.ax.set_ylabel('Density')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(4.5,7)
    plt.ylim(4.5,7)
    return ax



def project_scores(scores,seq):
    projected=scores*seq
    seq_length=seq.shape[0]
    absmax=np.argmax(np.absolute(projected),axis=1)
    projected=[projected[i][absmax[i]] for i in range(seq_length)]
    return projected

def shap_to_prob(projected_shap):
    projected_shap=np.squeeze(projected_shap)
    return np.abs(projected_shap)/np.sum(projected_shap) 

def get_vals_from_gkm_line(seq):
    return np.asarray([[float(i) for i in i.split(',')] for i in seq.split(';')])


def format_gkm_scores(gkm_file,chrom=None):
    pos_to_scores={}
    data=pd.read_csv(gkm_file,header=None,sep='\t')
    for index,row in data.iterrows():
        pos=row[0].split(':')
        cur_chrom=pos[0]
        if chrom is not None:
            if cur_chrom!=chrom:
                continue
        bp=int(pos[1])
        vals=get_vals_from_gkm_line(row[2])
        pos_to_scores[tuple([cur_chrom,bp])]=vals
    return pos_to_scores

def format_binary_deepshap(npz_file,chrom=None):
    data=np.load(npz_file)
    pos_to_scores={}
    bed_entries=data['bed_entries']
    interp_scores=data['interp_scores']
    num_regions=bed_entries.shape[0]
    for i in range(num_regions):
        cur_chrom=bed_entries[i][0]
        if chrom is not None:
            if cur_chrom!=chrom:
                continue
        cur_key=tuple([bed_entries[i][0],int(bed_entries[i][1])])
        cur_val=interp_scores[i][0]
        pos_to_scores[cur_key]=cur_val
    return pos_to_scores
    

