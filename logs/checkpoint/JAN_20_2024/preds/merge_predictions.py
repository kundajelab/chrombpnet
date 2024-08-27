import argparse
import os
import numpy as np
import pandas as pd
import subprocess
from pandas.errors import EmptyDataError
import deepdish as dd
import json
import time
import bigwig_helper
import pyfaidx
import h5py
import one_hot

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
folds = ["fold_0", "fold_1", "fold_2", "fold_3", "fold_4"]
#data = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",header=None)
data = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",header=None)
ddtpe="DNASE" 
cell_types = ["HEPG2"]

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)


def get_seq(peaks_df, genome, width=2114):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = []
    peaks_used = []
    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        vals.append(sequence)
    return one_hot.dna_to_one_hot(vals)


def write_predictions_h5py(output_prefix, logits, logcts, coords):

    # open h5 file for writing predictions
    output_h5_fname = "{}_predictions.h5".format(output_prefix)
    h5_file = h5py.File(output_h5_fname, "w")
    # create groups
    coord_group = h5_file.create_group("coords")
    pred_group = h5_file.create_group("predictions")

    num_examples=len(coords)

    coords_chrom_dset =  [str(coords[i][0]) for i in range(num_examples)]
    coords_start_dset =  [int(coords[i][1]) for i in range(num_examples)]
    coords_end_dset =  [int(coords[i][2]) for i in range(num_examples)]

    dt = h5py.special_dtype(vlen=str)

    # create the "coords" group datasets
    coords_chrom_dset = coord_group.create_dataset(
        "coords_chrom", data=np.array(coords_chrom_dset, dtype=dt),
        dtype=dt, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_start_dset", data=coords_start_dset, dtype=int, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_end_dset", data=coords_end_dset, dtype=int, compression="gzip")

    # create the "predictions" group datasets
    profs_dset = pred_group.create_dataset(
        "logits",
        data=logits,
        dtype=float, compression="gzip")
    logcounts_dset = pred_group.create_dataset(
        "logcounts", data=logcts,
        dtype=float, compression="gzip")

    # close hdf5 file
    h5_file.close()

#def _clip_prob(prob):
#    eps = np.finfo('float64').eps
#    return np.clip(prob, eps, 1 - eps)

#def convert_logits_to_logits(in_logits):
#	in_prob = softmax(in_logits)
#	clipped_prob = _clip_prob(in_prob)	
	
def merge_h5s_and_get_bigwig(genome, chrom_sizes, output_prefix):

	for cell_type in cell_types:
		output_prefix = output_prefix+ddtpe+"/"+cell_type+"/merge_folds_new/predictions_all_jan_2024"
		ndata = data[data[1]==cell_type].reset_index()
		for i,r in ndata.iterrows():
			print(i,r[2])

			temp_r="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/"+ddtpe+"/"+cell_type
			ofile=os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds_preds_done.txt")
			if os.path.isfile(ofile):
				h5_path = os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds_w_bias_predictions.h5")
				h5_path_corr = os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds_wo_bias_predictions.h5")

			elif os.path.isfile(os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds_preds_done.txt")):
				h5_path  = os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds_w_bias_predictions.h5")       
				h5_path_corr  = os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds_wo_bias_predictions.h5")       

			else:
				assert(False)
			           

			f = h5py.File(h5_path, "r")
			assert(f['predictions']['logits'].shape[0] == f['coords']['coords_chrom'].shape[0])
			print(type(f['coords']['coords_chrom'][:].tolist()))

			
			if i==0:
				coords_chrom_dset = f['coords']['coords_chrom'][:].tolist()
				coords_start_dset = f['coords']['coords_start_dset'][:].tolist()
				coords_end_dset = f['coords']['coords_end_dset'][:].tolist()

				mids = ((np.array(coords_start_dset)+np.array(coords_end_dset))//2).tolist()
				avg_logits = f['predictions']['logits'][:]
				avg_logcounts = f['predictions']['logcounts'][:]
				print(avg_logcounts[0:10])
				print(avg_logcounts.shape)
				print(avg_logits.shape)
		
			else:
				assert(f['coords']['coords_chrom'][:].tolist() == coords_chrom_dset)
				assert(f['coords']['coords_start_dset'][:].tolist() == coords_start_dset)
				assert(f['coords']['coords_end_dset'][:].tolist() == coords_end_dset)
				avg_logits += f['predictions']['logits'][:]
				avg_logcounts += f['predictions']['logcounts'][:]
				print(avg_logcounts[0:10])
				print(avg_logcounts.shape)
				print(avg_logits.shape)
		
		avg_logits = avg_logits / 5
		avg_logcounts = avg_logcounts / 5
		
		gs = bigwig_helper.read_chrom_sizes(chrom_sizes)
		
		regions=list(map(list,zip(coords_chrom_dset, coords_start_dset, coords_end_dset, mids)))
		write_predictions_h5py(output_prefix+"_w_bias", avg_logits, avg_logcounts, regions)
		final_array=softmax(avg_logits) * (np.expand_dims((np.exp(avg_logcounts)-1)[:,0],axis=1))
		print(final_array.shape)
		bigwig_helper.write_bigwig(final_array, 
							   regions, 
							   gs, 
							   output_prefix+"_w_bias.bw", 
							   outstats_file=output_prefix+"_w_bias.stat", 
							   debug_chr=False, 
							   use_tqdm=True)
							   
							   
		for i,r in ndata.iterrows():
			print(i,r[2])

			temp_r="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/"+ddtpe+"/"+cell_type
			ofile=os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds_preds_done.txt")
			if os.path.isfile(ofile):
				h5_path = os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds_w_bias_predictions.h5")
				h5_path_corr = os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds_wo_bias_predictions.h5")

			elif os.path.isfile(os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds_preds_done.txt")):
				h5_path  = os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds_w_bias_predictions.h5")       
				h5_path_corr  = os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds_wo_bias_predictions.h5")       

			else:
				assert(False)
			           
			f = h5py.File(h5_path_corr, "r")
			assert(f['predictions']['logits'].shape[0] == f['coords']['coords_chrom'].shape[0])
			print(type(f['coords']['coords_chrom'][:].tolist()))
			print(type(regions[0]))
			print(regions[0][0:10])

			assert(f['coords']['coords_chrom'][:].tolist() == coords_chrom_dset)
			assert(f['coords']['coords_start_dset'][:].tolist() == coords_start_dset)
			assert(f['coords']['coords_end_dset'][:].tolist() == coords_end_dset)

			if i==0:
				avg_logits = f['predictions']['logits'][:]
				avg_logcounts = f['predictions']['logcounts'][:]
				print(avg_logcounts[0:10])
				print(avg_logcounts.shape)
	
			else:

				avg_logits += f['predictions']['logits'][:]
				avg_logcounts += f['predictions']['logcounts'][:]
				print(avg_logcounts[0:10])
				print(avg_logcounts.shape)
		
		avg_logits = avg_logits / 5
		avg_logcounts = avg_logcounts / 5

		regions=list(map(list,zip(coords_chrom_dset, coords_start_dset, coords_end_dset, mids)))
		write_predictions_h5py(output_prefix+"_wo_bias", avg_logits, avg_logcounts, regions)
		final_array=softmax(avg_logits) * (np.expand_dims((np.exp(avg_logcounts)-1)[:,0],axis=1))
		print(final_array.shape)
		bigwig_helper.write_bigwig(final_array, 
							   regions, 
							   gs, 
							   output_prefix+"_wo_bias.bw", 
							   outstats_file=output_prefix+"_wo_bias.stat", 
							   debug_chr=False, 
							   use_tqdm=True)


		
if __name__=="__main__":


	genome="/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa"
	chrom_sizes="/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes"
	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/"
	merge_h5s_and_get_bigwig(genome, chrom_sizes, output_prefix)

	

