#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


dtype = "counts"
moods_concat="/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/"
#moods_concat="/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/"

moods_concat=moods_concat + "moods/"+  dtype + "/all_metaclusters.all_patterns.hits.csv"

moods_dict = pd.read_csv(moods_concat, header=None)

display(moods_dict.head())
print()
display(moods_dict.shape)

moods_dict['peak_chrom'] = moods_dict[0].apply(lambda x: x.split(':')[0])
moods_dict['peak_start'] = moods_dict[0].apply(lambda x: x.split(':')[1].split("-")[0]).astype(int)
moods_dict['peak_end'] = moods_dict[0].apply(lambda x: x.split(':')[1].split("-")[1]).astype(int)
display(moods_dict.head())


moods_dict['pattern_chr'] = moods_dict['peak_chrom']
moods_dict['pattern_start'] = moods_dict['peak_start'] + moods_dict[2].astype(int)
moods_dict['pattern_end'] = moods_dict['pattern_start'] + moods_dict[5].apply(lambda x: len(x)).astype(int)
moods_dict['pattern_name'] = moods_dict[1].apply(lambda x: x.split('.pfm')[0])

moods_dict = moods_dict[['pattern_chr', 'pattern_start', 'pattern_end', 'pattern_name', 4, 3, 5, 'peak_chrom', 'peak_start', 'peak_end']]
moods_dict.columns = ['pattern_chr', 'pattern_start', 'pattern_end', 'pattern_name', 'score', 'strand', 'sequence', 'peak_chrom', 'peak_start', 'peak_end']

display(moods_dict.head())
print()
display(moods_dict.shape)


# In[3]:


#tomtom = pd.read_csv("/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/SIGNAL/modisco_crop_500/" + dtype + ".tomtom.tsv", sep="\t")
tomtom = pd.read_csv("/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/HEPG2/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/SIGNAL/modisco_crop_500/" + dtype + ".tomtom.tsv", sep="\t")


# In[4]:


tomtom.head()


# In[5]:


label_dict = {}

for index,row in tomtom.iterrows():
    label_dict[row['Pattern']] = row['Match_1']
        

moods_dict['tf'] = moods_dict['pattern_name'].apply(lambda x: label_dict[x] if x in label_dict else 'tbd')
moods_dict = moods_dict[["pattern_chr", "pattern_start", "pattern_end", "tf", "pattern_name",
                         "score", "strand", "sequence", "peak_chrom", "peak_start", "peak_end"]]

display(moods_dict.head())
print()
display(moods_dict.shape)


# In[6]:


#main_dir = "/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/moods/" + dtype
main_dir = "/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/moods/" + dtype

moods_dict.drop_duplicates().to_csv(main_dir+"/moods_scan_on_modisco_full.bed", sep="\t", header=False, index=False)


# In[ ]:




