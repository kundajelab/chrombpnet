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
def make_plot(coord,
              label_prob,
              pred_prob,
              profile_shap,
              count_shap,
              label_sum,
              pred_sum,
              seq,
              binary_class_shap,
              binary_reg_shap,
              gkm_score,
              ymin=None,
              ymax=None,
              xmin=0,
              xmax=1000,
              motif_coords=[]):
    plt.rcParams["figure.figsize"]=25,10
    fig, axes = plt.subplots(6, 1)
    
    #predictions & labels from bpnet 
    axes[0].plot(label_prob,label='Label Prob',color='b')
    axes[0].plot(pred_prob,label='Pred Prob',color='r')
    axes[0].set_title(str(coord)+"; Count Label="+str(label_sum)+"; Count Pred="+str(pred_sum))   
    axes[0].legend() 
    axes[0].set_xlim(xmin,xmax)
    axes[0].set_xticks(list(range(xmin, xmax, 50,)))    
    
    #profile shap 
    axes[1]=plot_seq_importance(profile_shap,seq,xlim=(xmin,xmax),axes=axes[1])
    axes[1].set_title("BPNET Profile Loss SHAP")        
    if ymin is not None:
        axes[1].set_ylim(ymin,ymax)
    axes[1].set_xticks(list(range(xmin, xmax, 50,)))    

    #count shap 
    axes[2]=plot_seq_importance(count_shap,seq,xlim=(xmin,xmax),axes=axes[2])
    axes[2].set_title("BPNET Count Loss SHAP")
    if ymin is not None:
        axes[2].set_ylim(ymin,ymax)
    axes[2].set_xticks(list(range(xmin, xmax, 50,)))   

    #binary class shap 
    axes[3]=plot_seq_importance(binary_class_shap,seq,xlim=(xmin,xmax),axes=axes[3])
    axes[3].set_title("Binary Class. SHAP")
    if ymin is not None:
        axes[3].set_ylim(ymin,ymax)
    axes[3].set_xticks(list(range(xmin, xmax, 50,)))
    
    #binary regression shap 
    axes[4]=plot_seq_importance(binary_reg_shap,seq,xlim=(xmin,xmax),axes=axes[4])
    axes[4].set_title("Binary Reg. SHAP")
    if ymin is not None:
        axes[4].set_ylim(ymin,ymax)
    axes[4].set_xticks(list(range(xmin, xmax, 50,)))
    
    #GKM score 
    axes[5]=plot_seq_importance(gkm_score,seq,xlim=(xmin,xmax),axes=axes[5])
    axes[5].set_title("GKM Explain")
    if ymin is not None:
        axes[5].set_ylim(ymin,ymax)
    axes[5].set_xticks(list(range(xmin, xmax, 50,)))
    
    #add motif positions from vierstra 
    for motif_coord in motif_coords:
        for axes_index in range(6):
            axes[axes_index].axvspan(motif_coord[0], motif_coord[1], facecolor='b', alpha=0.3)
    plt.subplots_adjust(hspace=0.6)
    plt.show()

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
    

def get_motif_offsets_for_peaks(motif_intersections_path,chrom=None):
    '''
    motif_intersections_path is result of bedtools intersect on idr peak calls in our data and motif annotations (i.e. from vierstra) 
    '''
    motif_intersections=pd.read_csv(motif_intersections_path,header=None,sep='\t')    
    peak_tuple_to_motifs={}
    
    for index,row in motif_intersections.iterrows(): 
        peak_chrom=row[3]
        if chrom is not None: 
            if peak_chrom!=chrom:
                continue
        peak_start=row[4] 
        peak_summit_offset=row[6] 
        peak_summit=peak_start+peak_summit_offset
        region_start=peak_summit - 500
        region_end=peak_summit + 500
        #motif position 
        motif_offset_start=row[1]-region_start 
        if motif_offset_start<0: 
            continue 
        motif_offset_end=row[2]-region_start 
        if motif_offset_end > 1000:
            continue 
        cur_coord=tuple([peak_chrom,peak_summit])
        if cur_coord not in peak_tuple_to_motifs: 
            peak_tuple_to_motifs[cur_coord]=[[motif_offset_start,motif_offset_end]]
        else: 
            peak_tuple_to_motifs[cur_coord].append([motif_offset_start,motif_offset_end])
    return peak_tuple_to_motifs

