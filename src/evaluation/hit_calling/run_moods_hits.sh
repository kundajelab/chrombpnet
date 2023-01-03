# TFM_TFM_PATH = modisco hdf5 path
# TFM_TOMTOM_PATH = tomtom tsv path
# TFM_SHAP_PATH = interpretations bigwig
# TFM_PEAKS_PATH = regions for hit calling
# TFM_MOODS_DIR = output dir
# METHOD_INPUT = how to calculate contributions score for the hits? Look at Modes supported note below
# TFM_REFERENCE_PATH = reference path

### Modes supported METHOD_INPUT
### (1) sum = Take the sum of the contribution scores within hit window (preferes longer motifs)
### (2) mean = Take the mean of the contribution scores within hit window (prefers shorter motifs)
### (3) sum_norm = Take the sum of the contribution scores within hit window (preferes longer motifs) and divide by the mean of the total score in peak region
### (4) mean_norm = Take the mean of the contribution scores within hit window (prefers shorter motifs) and divide by the mean of the total score in peak region
### (3) and (4) Modes account for background - so you should get hits in regions with low overall importance
### I use mean_norm mode

### Make sure your peak regions (input to the argument TFM_PEAKS_PATH) are merged 
### Final output - overlaps_resolved_based_on_cwm_activations_normed.bed
### NOTE: The code outputs lot of intermediate files at every stage


export TFM_TFM_PATH=$1 \
export TFM_TOMTOM_PATH=$2 \
export TFM_SHAP_PATH=$3 \
export TFM_PEAKS_PATH=$4 \
export TFM_MOODS_DIR=$5 \
export METHOD_INPUT=$6 \
export TFM_REFERENCE_PATH=$7


jupyter nbconvert \
	--execute hit_calling.ipynb \
	--to HTML --output $TFM_MOODS_DIR/tfm_results \
	--ExecutePreprocessor.timeout=-1 &

wait

python break_overlaps_with_cwm_score.py -o $TFM_MOODS_DIR  -tfm $TFM_TFM_PATH -bw $TFM_SHAP_PATH -g $TFM_REFERENCE_PATH

wait

# you can ignore the below step if you dont want reports about co-occurence on the resolved overlaps hits file

jupyter nbconvert \
	--execute summarize_plots_on_resolved_overlaps.ipynb \
	--to HTML --output $TFM_MOODS_DIR/tfm_results_on_resolved_overlaps \
	--ExecutePreprocessor.timeout=-1 &

