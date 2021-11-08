import os
import math
import h5py
import modisco
import pandas as pd
import numpy as np
from IPython.core.display import HTML
from modisco.visualization import viz_sequence
from matplotlib import pyplot as plt
from plotnine import *

pd.options.display.max_colwidth = 500

setting="BIAS"
in_dataset=["K562"]
in_mode=["ATAC_10.14.2021_withinvivobias", "ATAC_10.14.2021_withinvivobias_500filts_mincount"]

input_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/'+setting
modisco_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/'+setting
tomtom_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_500/'+setting
viz_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper/visualize_tomtom/flank_len_500/'+setting
logo_link = 'http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper/modisco_logos/flank_len_500/'+setting
logo_dir_vier = '/oak/stanford/groups/akundaje/projects/chromatin-atlas/vierstra_logos/'
logo_link_vier = 'http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper/vierstra_logos/'

def path_to_image_html(path):
    return '<img src="'+ path + '" width="240" >'

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
    
def make_logo(match):
    if match + '.png' in os.listdir(logo_dir_vier):
        pass
    elif match + '.pfm' in os.listdir('/oak/stanford/groups/akundaje/soumyak/motifs/pfms'):
        ppm = np.loadtxt('/oak/stanford/groups/akundaje/soumyak/motifs/pfms/' + match + '.pfm', delimiter='\t')
        ppm = np.transpose(ppm)
        _plot_weights(viz_sequence.ic_scale(ppm, background=background),
                        path=logo_dir_vier + '/' + match + '.png')
        
def motif_dist(m1, cluster, tf, pattern):
    
    corrs = []
    
    m1 = np.loadtxt('/oak/stanford/groups/akundaje/soumyak/motifs/pfms/' + m1 + '.pfm', delimiter='\t')
    m1 = np.transpose(m1)
    m1_shape = m1.shape
    
    for x in ['fwd', 'rev']:
    
        m2 = np.genfromtxt(modisco_dir + '/' + cluster + '/' + tf + '/pattern_' + str(pattern) + '.' + x + '.ppm', delimiter='\t')
    
        m2_shape = m2.shape
        if m1_shape[0] > m2_shape[0]:
            diff = m1_shape[0] - m2_shape[0]
            for i in range(diff+1):
                new_m1 = m1[i:i+m2_shape[0]]
                dist = np.corrcoef(new_m1.ravel(), m2.ravel())[0,1]
                corrs.append(dist)
        else:
            diff = m2_shape[0] - m1_shape[0]
            for i in range(diff+1):
                new_m2 = m2[i:i+m1_shape[0]]
                dist = np.corrcoef(new_m2.ravel(), m1.ravel())[0,1]
                corrs.append(dist)
                
    max_corr = max(corrs)
    return max_corr

for dataset in in_dataset:
    if not os.path.isdir(viz_dir + '/' + dataset):
        os.mkdir(viz_dir + '/' + dataset)

    for mode in in_mode:
        print(mode)
        #mode = mode.replace("_new","")
        if not os.path.isdir(viz_dir + '/' + dataset + '/' + mode):
            os.mkdir(viz_dir + '/' + dataset + '/' + mode)
        for infile in os.listdir(input_dir + '/' + dataset + '/' + mode):
            print(infile)
            for score_type in ['count_shap', 'profile_shap']:
                print(score_type)
                print(tomtom_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.tsv')
                if os.path.isfile(tomtom_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.tsv'):
                    tomtom_file = tomtom_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.tsv'
                    print(tomtom_file)

                    tomtom_df = pd.read_csv(tomtom_file, sep='\t')
                    tomtom_df['modisco_cwm_fwd'] = [logo_link + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.pattern_' + str(i) + '.cwm.fwd.png'
                                                                      for i in range(len(tomtom_df))
                                                                     ]
                    tomtom_df['modisco_cwm_rev'] = [logo_link + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.pattern_' + str(i) + '.cwm.rev.png'
                                                                      for i in range(len(tomtom_df))
                                                                     ]

                    logo_dict = {x: [] for x in range(1,11)}

                    for index, row in tomtom_df.iterrows():
                        for i in range(1,11):
                            if not pd.isnull(row['match_' + str(i)]):
                                make_logo(row['match_' + str(i)])
                                logo_dict[i].append(logo_link_vier + '/' + row['match_' + str(i)] + '.png')
                            else:
                                logo_dict[i].append('NA')

                    for i in range(1,11):
                        tomtom_df['match' + str(i) + '_logo'] = logo_dict[i]

                    tomtom_df.columns = ['pattern', 'num_seqlets',
                                            'match0', 'qval0', 'match1', 'qval1',
                                            'match2', 'qval2', 'match3', 'qval3',
                                            'match4', 'qval4', 'match5', 'qval5',
                                            'match6', 'qval6', 'match7', 'qval7',
                                            'match8', 'qval8', 'match9', 'qval9',
                                            'modisco_cwm_fwd', 'modisco_cwm_rev',
                                            'match0_logo', 'match1_logo', 'match2_logo', 'match3_logo',
                                            'match4_logo', 'match5_logo', 'match6_logo', 'match7_logo',
                                            'match8_logo', 'match9_logo']

                    tomtom_df = tomtom_df[['pattern',
                                             'num_seqlets', 'modisco_cwm_fwd', 'modisco_cwm_rev',
                                             'match0', 'qval0', 'match0_logo', 'match1', 'qval1', 'match1_logo',
                                             'match2', 'qval2', 'match2_logo', 'match3', 'qval3', 'match3_logo',
                                             'match4', 'qval4', 'match4_logo', 'match5', 'qval5', 'match5_logo',
                                             'match6', 'qval6', 'match6_logo', 'match7', 'qval7', 'match7_logo',
                                             'match8', 'qval8', 'match8_logo', 'match9', 'qval9', 'match9_logo',
                                            ]]

                    tomtom_df.to_html(open(viz_dir + '/' + dataset + '/' + mode + '/' + infile + '.' + score_type + '.motifs.html', 'w'),
                              escape=False, formatters=dict(modisco_cwm_fwd=path_to_image_html,
                                                            modisco_cwm_rev=path_to_image_html,
                                                            match0_logo=path_to_image_html,
                                                            match1_logo=path_to_image_html,
                                                            match2_logo=path_to_image_html,
                                                            match3_logo=path_to_image_html,
                                                            match4_logo=path_to_image_html,
                                                            match5_logo=path_to_image_html,
                                                            match6_logo=path_to_image_html,
                                                            match7_logo=path_to_image_html,
                                                            match8_logo=path_to_image_html,
                                                            match9_logo=path_to_image_html
                                                            ), index=False)


tomtom_df.head()



