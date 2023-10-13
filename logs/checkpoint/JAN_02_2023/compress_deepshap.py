import deepdish as dd
import h5py 
import hdf5plugin
import numpy as np
import argparse
import copy
import subprocess
import os

def parse_args():
	parser = argparse.ArgumentParser(description="Compresses h5 file with SHAP scores and replaces the original file")
	parser.add_argument("-i", "--input_file", type=str, required=True, help="Input shap scores file to compress")
	parser.add_argument("-o", "--output_prefix", type=str, required=True, help="File name of the compressed output")	
	args = parser.parse_args()
	return args

def get_h5_type(h5_dict):
	return 0 if "shap" in h5_dict.keys() else 1


def compress_scores(h5_dict, h5_type, output_prefix):
	'''
	This function will take in a dictionary corresponding to the relevant h5 object.
	It will compress every array (as defined by the "seq" key) to np.float16
	'''
	print("Compressing")
	#For the chromatin atlas h5 type, we can just change the data type of the array
	if h5_type == 0:
		h5_compressed = copy.deepcopy(h5_dict)
		for key in h5_dict:
			data_vals = h5_dict[key]["seq"][:]
			if key != "raw":
				h5_compressed[key]["seq"] = data_vals.astype(np.float16)
			else:
				h5_compressed[key]["seq"] = data_vals.astype(np.int8)
		dd.io.save(output_prefix + "_compressed.h5", h5_compressed, compression="blosc")
	#For the tf-atlas h5 type, we create a new dataset for each array with the appropriate type
	elif h5_type == 1:
		h5_compressed = h5py.File(output_prefix + "_compressed.h5", "a")
		num_examples = h5_dict["hyp_scores"].shape[0]
		seq_len = 2114		
		coords_chrom_dset = h5_compressed.create_dataset(
		    "coords_chrom", (num_examples,),
		    dtype=h5py.string_dtype(encoding="ascii"), 
		    **hdf5plugin.Blosc()
		)
		coords_chrom_dset[:] = h5_dict['coords_chrom'][:].astype('U8')

		coords_start_dset = h5_compressed.create_dataset(
		    "coords_start", (num_examples,), dtype="i4", 
		    **hdf5plugin.Blosc()
		)
		coords_start_dset[:] = h5_dict['coords_start'][:]

		coords_end_dset = h5_compressed.create_dataset(
		    "coords_end", (num_examples,), dtype="i4", 
		    **hdf5plugin.Blosc()
		)
		coords_end_dset[:] = h5_dict['coords_end'][:]

		hyp_scores_dset = h5_compressed.create_dataset(
		    "hyp_scores", (num_examples, seq_len, 4), dtype="f2",
		    **hdf5plugin.Blosc()
		)
		hyp_scores_dset[:, :, :] = h5_dict['hyp_scores'][:].astype(np.float16)

		input_seqs_dset = h5_compressed.create_dataset(
		    "input_seqs", (num_examples, seq_len, 4), dtype="i1",
		    **hdf5plugin.Blosc()
		)
		input_seqs_dset[:, :, :] = h5_dict['input_seqs'][:].astype(np.int8)
		h5_compressed.close()

def quality_check(h5_compressed, h5_orig, h5_type):
	'''
	Here, we ensure the compressed h5 data has been produced without any errors by comparing it to the original data
	'''
	print("Running quality check")
	#We flatten the SHAP arrays and check the length and the correlations look good
	if h5_type == 0:
		comp_flat = h5_compressed["shap"]["seq"].flatten()
		orig_flat = h5_orig["shap"]["seq"].flatten()
	elif h5_type == 1:
		comp_flat = h5_compressed["hyp_scores"][:].flatten()
		orig_flat = h5_orig["hyp_scores"][:].flatten()	

	assert len(comp_flat) == len(orig_flat), "Compressed file does not have the full set of values"
	assert np.corrcoef(comp_flat, orig_flat)[0,1] > 0.98, "Compressed file is not sufficiently similar to the original file"

def modisco_check(h5_compressed, output_prefix, h5_type):
	'''
	After we have reloaded the saved compressed file, we do a test modisco run just to make sure there are no issues
	The run is not supposed to produce anything informative; we just want to make sure it runs without any errors. 
	'''
	print("Running test modisco run")
	#Selected_samples is the first 5 shap scores
	#seq_array is the corresponding DNA sequences
	if h5_type == 0:
		selected_samples = h5_compressed["shap"]["seq"][:5]
		if "raw" in h5_compressed.keys():
			seq_array = h5_compressed["raw"]["seq"][:5]
		else:
			seq_array = np.eye(4)[np.random.choice([0,1,2,3], [selected_samples.shape[0], selected_samples.shape[2]])].transpose(0,2,1).astype(np.int8)
	elif h5_type == 1:
		selected_samples = h5_compressed["hyp_scores"][:][:5].transpose(0,2,1)
		if "input_seqs" in h5_compressed.keys():
			seq_array = h5_compressed["input_seqs"][:][:5].transpose(0,2,1)
		else:
			seq_array = np.eye(4)[np.random.choice([0,1,2,3], [selected_samples.shape[0], selected_samples.shape[2]])].transpose(0,2,1).astype(np.int8)
	#We save the two arrays as npy files and run modisco
	print(selected_samples.shape, seq_array.shape)
	np.save(output_prefix+"_test_modisco_shap.npy", selected_samples)
	np.save(output_prefix+"_test_modisco_seqs.npy", seq_array)
	modisco_command = ["modisco", "motifs", "-s", output_prefix+"_test_modisco_seqs.npy", "-a", output_prefix+"_test_modisco_shap.npy", "-n", "200", "-o", output_prefix+"_modisco_results.h5"]
	modisco_run = subprocess.run(modisco_command)
	assert modisco_run.returncode == 0, "Running modisco returned an error"
	os.remove(output_prefix+"_test_modisco_shap.npy")
	os.remove(output_prefix+"_test_modisco_seqs.npy")
	os.remove(output_prefix+"_modisco_results.h5")

def main():
	args = parse_args()
	#We first test what fields are in the h5 file. If it is the type with "shap" as one of the fields, we can open with deepdish, and it's easier
	#If not, we open with h5py instead
	h5_type_test = h5py.File(args.input_file, "r")
	h5_type = get_h5_type(h5_type_test)
	h5_type_test.close()
	#The procedure works as follows:
	#1. We save a new h5 file with the appropriate compression
	#2. We reload the compressed file and perform quality checks on it
	#3. We run modisco on the compressed data and ensure it works without errors
	if h5_type == 0:
		h5_orig = dd.io.load(args.input_file)
		compress_scores(h5_orig, h5_type, args.output_prefix)
		h5_compressed_reloaded = dd.io.load(args.output_prefix + "_compressed.h5")
	else:
		h5_orig = h5py.File(args.input_file, "r")
		compress_scores(h5_orig, h5_type, args.output_prefix)
		h5_compressed_reloaded = h5py.File(args.output_prefix + "_compressed.h5", "r")
	quality_check(h5_compressed_reloaded, h5_orig, h5_type)
	modisco_check(h5_compressed_reloaded, args.output_prefix, h5_type)
	#os.remove(args.input_file)




if __name__ == "__main__":
	main()
