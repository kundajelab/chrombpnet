import re

import h5py
import hdf5plugin
import numpy as np

import argparse
import gc

from typing import List, Optional
from pathlib import Path

import logging

import pandas as pd
import matplotlib.pyplot as plt

import psutil


def plot_shaps(
        shap_list: List[np.ndarray], n_points: int, plot_title: Optional[str], output_dir: Path
) -> None:
    """Takes in mean_shap values and plots a random number of points"""

    rand_generator = np.random.default_rng(seed=123)
    max_points = shap_list[0].shape[0] * shap_list[0].shape[1] * shap_list[0].shape[2]
    if n_points < max_points:
        sampling_idx = rand_generator.integers(low=0, high=max_points, size=n_points)
    else:
        sampling_idx = range(0, max_points + 1)

    samples_dict = {}
    for idx, shap in enumerate(shap_list):
        # flatten shaps, then sample them
        fold = f"fold_{idx}"
        samples_dict[fold] = shap.flatten()[sampling_idx]

    # convert to data frame
    shap_df = pd.DataFrame(data=samples_dict)

    # plot shap_df in pairwise scatter
    g = pd.plotting.scatter_matrix(shap_df)
    if plot_title is not None:
        plt.suptitle(plot_title)

    plt.savefig(output_dir / f"interpret_{plot_title}.png")


def average_across_shap_files(
        shaps_list: List[Path],
        output_path: Path,
        plot_data: bool,
        plot_title: Optional[str],
        n_points: int
):
    logging.basicConfig(level=logging.NOTSET)
    logger = logging.getLogger(__name__)

    # check that length of file list is as expected
    assert len(shaps_list) == 5, f"Number of files to average is {len(shaps_list)}"

    # iterate over files in shaps_list, one key at a time
    for key in ['coords_start', 'coords_chrom', 'coords_end', 'raw', 'shap', 'projected_shap']:

        # initialize one list to store the values in the h5 file, one key at a time
        key_list: List[np.array] = []

        for shap_h5 in shaps_list:
            assert shap_h5.is_file(), f"{shap_h5} does not exist"
            logger.info(f'Loading shap files for {shap_h5} and {key}...')
            with h5py.File(shap_h5, 'r') as f_in:
                assert key in f_in.keys(), f"{key} does not exist!"
                if key == 'coords_chrom':
                    key_list.append([c.decode("utf-8") for c in f_in["coords_chrom"][:]])
                elif 'coords' in key:
                    key_list.append(f_in[key][()])
                else:
                    key_list.append(f_in[key]['seq'][()])

        # if it's the raw key, need to make sure that all seqs are the same, then
        # write just the first element of the list to file
        # note that here we use 'w' rather than 'a' -- this is in case a file already exists in
        # output_path (we will write over it)
        # add coords keys so that don't have to do this after the fact
        if key == 'coords_start':
            compare_seqs = [np.array_equal(a, b) for a in key_list for b in key_list]
            assert all(compare_seqs), "Not all start values are equal! Please check regions inputs to shap calculations"
            with h5py.File(output_path, 'w') as f:
                logger.info(f'Saving {key} to file...')
                f.create_dataset(name=f'{key}', data=key_list[0], **hdf5plugin.Blosc())

        elif key == 'coords_end':
            compare_seqs = [np.array_equal(a, b) for a in key_list for b in key_list]
            assert all(compare_seqs), "Not all end values are equal! Please check regions inputs to shap calculations"
            with h5py.File(output_path, 'a') as f:
                logger.info(f'Saving {key} to file...')
                f.create_dataset(name=f'{key}', data=key_list[0], **hdf5plugin.Blosc())

        elif key == 'coords_chrom':
            compare_seqs = [np.array_equal(a, b) for a in key_list for b in key_list]
            assert all(compare_seqs), "Not all chrom values are equal! Please check regions inputs to shap calculations"
            with h5py.File(output_path, 'a') as f:
                logger.info(f'Saving {key} to file...')
                f.create_dataset(name=f'{key}', data=key_list[0])

        elif key == 'raw':
            with h5py.File(output_path, 'a') as f:
                logger.info(f'Saving {key} to file...')
                f.create_dataset(name=f'{key}/seq', data=key_list[0], **hdf5plugin.Blosc())

        # if it's the projected_shap key, decide whether to plot the data, and then save the mean across shaps
        elif key == 'projected_shap':
            if plot_data:
                logger.info('Plotting data...')
                plot_shaps(shap_list=key_list, n_points=n_points, plot_title=plot_title, output_dir=output_path.parent)
            with h5py.File(output_path, 'a') as f:
                logger.info(f'Saving {key} to file...')
                logger.info(f"Memory used: {psutil.virtual_memory()[3] / 10 ** 9}")
                f.create_dataset(
                    name=f'{key}/seq',
                    data=np.nanmean(np.array(key_list), axis=0, dtype=np.float16),
                    **hdf5plugin.Blosc()
                )

        # otherwise, just take the mean across values and save them
        elif key == 'shap':
            with h5py.File(output_path, 'a') as f:
                logger.info(f'Saving {key} to file...')
                logger.info(f"Memory used: {psutil.virtual_memory()[3] / 10 ** 9}")
                f.create_dataset(
                    name=f'{key}/seq',
                    data=np.nanmean(np.array(key_list), axis=0, dtype=np.float16),
                    **hdf5plugin.Blosc()
                )
        # not sure if this helps, but might be worth a shot
        del key_list
        gc.collect()
        mem_used = psutil.virtual_memory()[3]
        logger.info(f"Memory used: {mem_used / 10 ** 9}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="calculate mean shap over the given h5s")
    parser.add_argument("--in_dir", type=str,
                        required=True, help="Directory(ies) with counts and profile shaps", nargs='+')
    parser.add_argument("--output_dir", type=str, default='/cromwell_root/',
                        help="output directory for the mean_shap")
    parser.add_argument("--plot_data", type=bool, default=False)
    parser.add_argument("--n_points", type=int, default=20000)
    parser.add_argument("--experiment", type=str, default="Experiment", help="Annotation for plot title (e.g. encid)")

    args = parser.parse_args()
    # also make sure to exclude shaps that were calculated as part of QC subsampling
    counts_shaps = []
    profile_shaps = []
    # shaps may be spread across multiple directories; find all of them
    for cur_path in args.in_dir:
        cur_counts_shaps = [x for x in list(Path(cur_path).glob('**/*.counts_scores.h5'))
                            if not (re.search("subsample", str(x)))]
        counts_shaps += cur_counts_shaps
        cur_profile_shaps = [x for x in list(Path(cur_path).glob('**/*.profile_scores.h5'))
                             if not (re.search("subsample", str(x)))]
        profile_shaps += cur_profile_shaps

    average_across_shap_files(
        shaps_list=counts_shaps,
        output_path=Path(f"{args.output_dir}/counts_mean_shap_scores.h5"),
        plot_data=args.plot_data,
        plot_title=f"{args.experiment}_counts_shaps",
        n_points=args.n_points
    )

    average_across_shap_files(
        shaps_list=profile_shaps,
        output_path=Path(f"{args.output_dir}/profile_mean_shap_scores.h5"),
        plot_data=args.plot_data,
        plot_title=f"{args.experiment}_profile_shaps",
        n_points=args.n_points
    )

