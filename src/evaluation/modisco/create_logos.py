
import os
import h5py
import modisco
import numpy as np
from modisco.visualization import viz_sequence
from matplotlib import pyplot as plt

setting="SIGNAL"
in_dataset=["K562"]
in_mode=["4_4_shifted_ATAC_09.12.2021_bias_filters_500_new"]

input_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/'+setting
modisco_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/'+setting
logo_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco_logos/flank_len_500/'+setting


def _plot_weights(array,
                  path,
                  figsize=(10,3),
                 **kwargs):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111) 
    viz_sequence.plot_weights_given_ax(ax=ax, array=array,**kwargs)
    plt.savefig(path)
    plt.close()


background = np.array([0.25, 0.25, 0.25, 0.25])
for dataset in in_dataset:
    if not os.path.isdir(logo_dir + '/' + dataset):
        os.mkdir(logo_dir + '/' + dataset)
    
    for mode in in_mode:
        print(mode)
        #mode = mode.replace("_new","")
        if not os.path.isdir(logo_dir + '/' + dataset + '/' + mode):
            os.mkdir(logo_dir + '/' + dataset + '/' + mode)
        for infile in os.listdir(input_dir + '/' + dataset + '/' + mode):
            print(infile)
            for score_type in ['count_shap', 'profile_shap']:
                print(score_type)
                if os.path.isfile(modisco_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.hdf5'):
                    modisco_file = modisco_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.hdf5'
                    print(modisco_file)
                    hdf5_results = h5py.File(modisco_file,'r')
                    for metacluster_name in hdf5_results["metacluster_idx_to_submetacluster_results"]:
                        metacluster = hdf5_results["metacluster_idx_to_submetacluster_results"][metacluster_name]
                        if metacluster['activity_pattern'][0] == 1:
                            all_pattern_names = [x.decode("utf-8") for x in list(metacluster["seqlets_to_patterns_result"]["patterns"]["all_pattern_names"][:])]
                            for pattern_name in all_pattern_names:
                                cwm_fwd = np.array(metacluster['seqlets_to_patterns_result']['patterns'][pattern_name]['task0_contrib_scores']['fwd'])
                                cwm_rev = np.array(metacluster['seqlets_to_patterns_result']['patterns'][pattern_name]['task0_contrib_scores']['rev'])

                                score_fwd = np.sum(np.abs(cwm_fwd), axis=1)
                                score_rev = np.sum(np.abs(cwm_rev), axis=1)

                                trim_thresh_fwd = np.max(score_fwd) * 0.3
                                trim_thresh_rev = np.max(score_rev) * 0.3

                                pass_inds_fwd = np.where(score_fwd >= trim_thresh_fwd)[0]
                                pass_inds_rev = np.where(score_rev >= trim_thresh_rev)[0]

                                start_fwd, end_fwd = max(np.min(pass_inds_fwd) - 4, 0), min(np.max(pass_inds_fwd) + 4 + 1, len(score_fwd) + 1)
                                start_rev, end_rev = max(np.min(pass_inds_rev) - 4, 0), min(np.max(pass_inds_rev) + 4 + 1, len(score_rev) + 1)

                                trimmed_cwm_fwd = cwm_fwd[start_fwd:end_fwd]
                                trimmed_cwm_rev = cwm_rev[start_rev:end_rev]

                                _plot_weights(trimmed_cwm_fwd,
                                             path=logo_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.' + pattern_name + '.cwm.fwd.png')
                                _plot_weights(trimmed_cwm_rev,
                                             path=logo_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.' + pattern_name + '.cwm.rev.png')


# In[ ]:




