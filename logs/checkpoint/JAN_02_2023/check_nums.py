import pandas as pd

predictions_bed = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/IMR90_new/preds_upload/fold_0/IMR90_new_w_bias_all_regions.bed", sep='\t', header=None).drop_duplicates()
predictions_bed = predictions_bed[predictions_bed[0] != "chrM"].reset_index(drop=True)

interpret_bd = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/IMR90_new/merge_folds_all_regions_may_05_24/IMR90_new_folds_merged.profile_scores_new_compressed.bed.gz", sep='\t', header=None)
interpret_bd[10] = interpret_bd[1]+interpret_bd[9]-500
interpret_bd[11] = interpret_bd[1]+interpret_bd[9]+500

interpret_bd = interpret_bd[[0,10,11]].drop_duplicates().rename(columns={10:1, 11:2}).sort_values(by=[0,1,2]).reset_index(drop=True)
predictions_bed = predictions_bed.sort_values(by=[0,1,2]).reset_index(drop=True)


print(predictions_bed.head())
print(interpret_bd.head())

print(predictions_bed.equals(interpret_bd))

import pybedtools

x = pybedtools.BedTool.from_dataframe(predictions_bed)
x = x.merge()
y = x.to_dataframe()

print(y.head())

temp="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/IMR90_new/preds_upload/average_preds/"
y.to_csv(temp+"merged.viz.bed.gz", compression='gzip', sep='\t', header=False, index=False)

print(interpret_bd.shape)
print(predictions_bed.shape)
print(len(set(interpret_bd[0])))
print(len(set(predictions_bed[0])))


