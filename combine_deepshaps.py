import argparse
import deepdish
import os
import pandas as pd
import deepdish as dd
import numpy as np
def parse_args():
        parser = argparse.ArgumentParser(description="Compresses h5 file with SHAP scores and replaces the original file")
        parser.add_argument("-i", "--input_dir", type=str, required=True, help="")
        parser.add_argument("-if", "--input_file", type=str, nargs='+', required=True, help="")
        parser.add_argument("-e", "--experiment", type=str, required=True, help="")
        parser.add_argument("-o", "--output_file", type=str, required=True, help="")
        parser.add_argument("-t", "--type", type=str, required=True, help="")
        args = parser.parse_args()

        return args

def main(args):

	interpretation_files = []
	bed_files = []
	for filer in args.input_file:
		interpret_f = os.path.join(args.input_dir,os.path.join(filer,"interpret/full_"+args.experiment+"."+args.type+"_scores.h5"))
		if os.path.isfile(interpret_f):
			interpretation_files.append(interpret_f)
		bed_file = os.path.join(args.input_dir,os.path.join(filer,"interpret/full_"+args.experiment+".interpreted_regions_"+args.type+".bed"))
		if os.path.isfile(bed_file):
			bed_files.append(bed_file)
	assert(len(interpretation_files)==5)
	assert(len(bed_files)==5)

	idx = 0
	for bed_file in bed_files:
		print(bed_file)
		if idx==0:
			bed1 = pd.read_csv(bed_file,sep="\t", header=None)
		else:
			bed  = pd.read_csv(bed_file,sep="\t", header=None)
			assert(bed1.equals(bed))
		idx+=1
	idx = 0
	for interpret_f in interpretation_files:
		print(interpret_f)
		data = dd.io.load(interpret_f)["shap"]["seq"]	
		if idx==0:
			new_d = data
		else:
			assert(new_d.shape==data.shape)
			assert(new_d.shape[0]==bed.shape[0])
			new_d = new_d + data
			#print(data[0][0][0])
		del data
		#print(new_d.shape)
		#print(bed.shape)
		#print(new_d[0][0][0])
		idx += 1

	new_d = new_d / 5
	#print(new_d[0][0][0])

	dd.io.save("{}.mean_shap.".format(args.output_file)+args.type+"_scores.h5", new_d, compression='blosc')

if __name__=="__main__":
	args = parse_args()
	main(args)	
