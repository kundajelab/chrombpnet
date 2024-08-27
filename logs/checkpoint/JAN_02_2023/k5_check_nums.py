import pandas as pd

predictions_bed = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/preds_upload/fold_0/K562_w_bias_all_regions.bed", sep='\t', header=None).drop_duplicates()
predictions_bed = predictions_bed[predictions_bed[0] != "chrM"].reset_index(drop=True)
predictions_bed = predictions_bed[~predictions_bed[0].str.contains("random")].reset_index(drop=True)
predictions_bed = predictions_bed[~predictions_bed[0].str.contains("EMV")].reset_index(drop=True)



#interpret_bd = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/merge_folds_new_may_05_24/K562_folds_merged.profile_scores_new_compressed.bed", sep='\t', header=None, compression='gzip')
interpret_bd = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/k562.merged.atac.dnase.peaks.bed", sep='\t', header=None)
interpret_bd[10] = interpret_bd[1]+interpret_bd[9]-500
interpret_bd[11] = interpret_bd[1]+interpret_bd[9]+500

interpret_bd = interpret_bd[~interpret_bd[0].str.contains("random")].reset_index(drop=True)
interpret_bd = interpret_bd[~interpret_bd[0].str.contains("EMV")].reset_index(drop=True)

interpret_bd = interpret_bd[[0,10,11]].drop_duplicates().rename(columns={10:1, 11:2}).sort_values(by=[0,1,2]).reset_index(drop=True)
predictions_bed = predictions_bed.sort_values(by=[0,1,2]).reset_index(drop=True)


print(predictions_bed.head())
print(interpret_bd.head())

print(predictions_bed.shape)
print(interpret_bd.shape)
print(len(interpret_bd.merge(predictions_bed).drop_duplicates()))
print(len(interpret_bd.drop_duplicates()))
#print(predictions_bed.equals(interpret_bd))
print(predictions_bed.merge(interpret_bd).drop_duplicates().equals(predictions_bed.drop_duplicates()))


import pybedtools

x = pybedtools.BedTool.from_dataframe(predictions_bed)
x = x.merge()
y = x.to_dataframe()

print(y.shape)

#x = pybedtools.BedTool.from_dataframe(interpret_bd)
#x = x.merge()
#y = x.to_dataframe()

#print(y.shape)


#temp="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/preds_upload/average_preds/"
#y.to_csv(temp+"merged.viz.bed.gz", compression='gzip', sep='\t', header=False, index=False)

print(interpret_bd.shape)
print(predictions_bed.shape)
print(len(set(interpret_bd[0])))
print(len(set(predictions_bed[0])))
print(set(interpret_bd[0]) == set(predictions_bed[0]))

