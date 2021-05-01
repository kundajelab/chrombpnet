import argparse
import pdb 
import pysam 
import numpy as np
import pandas as pd
from plot_letters import * 
from kerasAC.splits import *
import math
from math import floor
import matplotlib
from matplotlib import pyplot as plt 

ltrdict = {'a':[1,0,0,0],
           'c':[0,1,0,0],
           'g':[0,0,1,0],
           't':[0,0,0,1],
           'n':[0,0,0,0],
           'A':[1,0,0,0],
           'C':[0,1,0,0],
           'G':[0,0,1,0],
           'T':[0,0,0,1],
           'N':[0,0,0,0]}


def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--ref_fasta",default="/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta")
    parser.add_argument("--bed")
    parser.add_argument("--scores_prefix")
    parser.add_argument("--flank",type=int,default=500) 
    parser.add_argument("--outf_prefix")
    parser.add_argument("--plot_start_base",type=int,default=400)
    parser.add_argument("--plot_end_base",type=int,default=600)
    return parser.parse_args()

def plot_seq_importance(outf,tracks,labels,ylim,xlim, heatmap_indices=None, figsize=(20,20)):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 10}
    matplotlib.rc('font', **font)
    num_plots=len(tracks)
    f,axes=plt.subplots(nrows=10,ncols=1,dpi=100,figsize=figsize)
    show=False
    seq_len = tracks[0].shape[0]
    for plot_index in range(num_plots): 
        cur_track=tracks[plot_index]
        cur_ylim=ylim[plot_index]
        cur_xlim=xlim[plot_index]
        vmin=-1*max([abs(cur_ylim[0]),abs(cur_ylim[1])])
        vmax=max([abs(cur_ylim[0]),abs(cur_ylim[1])])
        axes[plot_index]=plot_bases_on_ax(cur_track,axes[plot_index],show_ticks=True)
        axes[plot_index].set_xlim(cur_xlim) 
        axes[plot_index].set_ylim(cur_ylim)
        cur_label=labels[plot_index] 
        axes[plot_index].set_title(cur_label)
        axes[plot_index].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
    plt.subplots_adjust(hspace=0.5)
    plt.tight_layout()
    plt.savefig(outf,format='png',dpi=120)
    plt.close() 
    return 

def get_vals_from_gkm_line(seq):
    return np.asarray([[float(i) for i in i.split(',')] for i in seq.split(';')])

def one_hot_encode(seq):
    return np.array([ltrdict.get(x,[0,0,0,0]) for x in seq])



def main():
    args=parse_args()
    #load the gkmexplain data
    regions=pd.read_csv(args.bed,header=0,sep='\t')
    regions_to_vals={}
    ref=pysam.FastaFile(args.ref_fasta)
    for fold in range(10):
        scores=pd.read_csv(args.scores_prefix+'.'+str(fold),header=None,sep='\t')
        print("loaded gkmexplain scores") 
        for index,row in scores.iterrows():
            cur_region_info=row[0]
            chrom=row[0].split(':')[0]
            start_pos=int(row[0].split(':')[1].split('-')[0])
            end_pos=int(row[0].split(':')[1].split('-')[1])
            pred=row[1]
            vals=get_vals_from_gkm_line(row[2])
            seq=one_hot_encode(ref.fetch(chrom,start_pos,end_pos))
            if cur_region_info not in regions_to_vals: 
                regions_to_vals[cur_region_info]={}
            if fold not in regions_to_vals[cur_region_info]:
                regions_to_vals[cur_region_info][fold]={} 
            regions_to_vals[cur_region_info][fold]['vals']=vals
            regions_to_vals[cur_region_info][fold]['pred']=pred
            regions_to_vals[cur_region_info][fold]['seq']=seq


    for region in regions_to_vals:
        plot_wrapper(region,
                     vals=[regions_to_vals[region][i]['vals'] for i in range(10)],
                     seq=[regions_to_vals[region][i]['seq'] for i in range(10)],
                     preds=[regions_to_vals[region][i]['pred'] for i in range(10)],
                     args=args)
        
    
def plot_wrapper(region,
                 vals,
                 seq,
                 preds,
                 args):
    png_title=args.outf_prefix+'.'+region+'.allfolds.png'
    toplot_tracks=[]
    toplot_labels=[]
    for i in range(10):
        track=vals[i]*seq[i]
        toplot_tracks.append(track)
        toplot_labels.append(region+' fold '+str(i)+' prediction:'+str(round(preds[i],2)))
        
    minvals=[]
    maxvals=[]
    
    #gkm y bounds 
    gkm_min=min([np.amin(i) for i in toplot_tracks])
    minvals=[gkm_min]*len(toplot_tracks)
    gkm_max=max([np.amax(i) for i in toplot_tracks])
    maxvals=[gkm_max]*len(toplot_tracks)

    ylim=[(minvals[i],maxvals[i]) for i in range(len(toplot_tracks))]
    xlim=[(args.plot_start_base,args.plot_end_base) for i in range(len(toplot_tracks))]
          
    plot_seq_importance(png_title,
                        toplot_tracks,
                        toplot_labels,
                        ylim=ylim,
                        xlim=xlim,
                        heatmap_indices=None)
    
if __name__=="__main__":
    main()
    
