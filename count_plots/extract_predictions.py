import  h5py
import pandas as pd
import numpy as np 
import argparse
import pdb

import matplotlib 
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

plt.rcParams["figure.figsize"]=5,4
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}

matplotlib.rc('font', **font)
#Performance metrics for profile models



def parse_args():
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("-predictions_files",nargs="+")
    parser.add_argument("-formats",nargs="+",choices=['legacy','new'])
    parser.add_argument("-labels",nargs="+") 
    parser.add_argument("-outf")
    return parser.parse_args() 

def parse_new_format(curfile):
    data=h5py.File(curfile,'r')
    coords=data['coords'][:]
    coords=[[i.decode('utf8')  for i in j] for j in coords]
    coords=[tuple([i[0],int(i[1])]) for i in coords]
    preds=data['pred_1'][:]
    cur_dict={}
    for i in range(len(coords)):
        cur_dict[tuple(coords[i][0:2])]=preds[i]
    return pd.DataFrame.from_dict(cur_dict,orient='index')


def parse_legacy_format(curfile):
    data=pd.DataFrame.from_dict(pd.read_hdf(curfile).to_dict()[0],orient='index')    
    return data 

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

    ax.scatter( x, y, c=z )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Density')
    plt.plot([i for i in range(14)],[i for i in range(14)],'k')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(0,14)
    plt.ylim(0,14)
    return ax



def main():
    args=parse_args()
    tomerge=[] 
    for i in range(len(args.predictions_files)):
        cur_f=args.predictions_files[i]
        cur_format=args.formats[i]
        if cur_format=="new":
            cur_dict=parse_new_format(cur_f)
        else:
            assert cur_format=="legacy"
            cur_dict=parse_legacy_format(cur_f)
        tomerge.append(cur_dict)
    merged=pd.merge(left=tomerge[0],right=tomerge[1],how='outer',left_index=True,right_index=True)
    merged.rename(columns={'0_x':args.labels[0],'0_y':args.labels[1]},inplace=True)
    merged.to_csv(args.outf+".tsv",header=True,index=True,sep='\t')
    print("wrote outputs")
    merged=merged.dropna()
    density_scatter(merged[args.labels[0]].values,
                    merged[args.labels[1]].values,
                    xlab=args.labels[0],
                    ylab=args.labels[1])
    plt.savefig(args.outf+'.png',format='png',dpi=300)

if __name__=="__main__":
    main()
    
