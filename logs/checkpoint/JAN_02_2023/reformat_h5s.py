import deepdish as dd
import argparse
import pandas as pd
import h5py
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-h5py1", "--h5_path1", type=str)
parser.add_argument("-r1", "--regions1", type=str)
parser.add_argument("-o", "--output_prefix", type=str)
args = parser.parse_args()



def write_predictions_h5py(output_prefix, scores, coords_chrom_dset, coords_start_dset, coords_end_dset):

    # open h5 file for writing predictions
    output_h5_fname = "{}_attribs_reformatted.h5".format(output_prefix)
    h5_file = h5py.File(output_h5_fname, "w")
    # create groups
    coord_group = h5_file.create_group("coords")
    pred_group = h5_file.create_group("attributions")

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
        "shap",
        data=scores,
        dtype=float, compression="gzip")

    # close hdf5 file
    h5_file.close()

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit", "temp"]

try:
    regions_df = pd.read_csv(args.regions1, sep="\t", header=None, names=NARROWPEAK_SCHEMA, compression='gzip')
except:
    regions_df = pd.read_csv(args.regions1, sep="\t", header=None, names=NARROWPEAK_SCHEMA)

regions_df["temp_start"] = regions_df["start"]+regions_df["summit"]-(2114//2)
regions_df["temp_end"] = regions_df["start"]+regions_df["summit"]+(2114//2)
regions_df["temp_sum"] = regions_df["start"]+regions_df["summit"]
regions = regions_df[["chr", "temp_start", "temp_end", "temp_sum"]].values.tolist()
print(regions_df.head())
num_examples = len(regions)
print(num_examples)

coords_chrom_dset1 =  [str(regions[i][0]) for i in range(num_examples)]
coords_start_dset1 =  [int(regions[i][1]) for i in range(num_examples)]
coords_end_dset1 =  [int(regions[i][2]) for i in range(num_examples)]

scores = dd.io.load(args.h5_path1)
print(scores['shap']['seq'].shape)
assert(scores['shap']['seq'].shape[0]==num_examples)

write_predictions_h5py(args.output_prefix, scores['shap']['seq'], coords_chrom_dset1, coords_start_dset1, coords_end_dset1)
