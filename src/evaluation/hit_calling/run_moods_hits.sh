dtype="counts"
#TFM_TFM_PATH = modisco hdf5 path
#TFM_TOMTOM_PATH = tomtom tsv path
#TFM_SHAP_PATH = interpretations bigwig
#TFM_PEAKS_PATH = regions for hit calling
#TFM_MOODS_DIR = output dir
#METHOD_INPUT = how to calculate contributions score? - 
#TFM_REFERENCE_PATH = reference path


export TFM_TFM_PATH=$1/modisco_results_allChroms_$dtype.hdf5 \
export TFM_TOMTOM_PATH=$1/$dtype.tomtom.tsv \
export TFM_SHAP_PATH=$2.counts.bw \
export TFM_PEAKS_PATH=$4 \
export TFM_MOODS_DIR=$3 \
export METHOD_INPUT=$5 \
export TFM_REFERENCE_PATH=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa 


jupyter nbconvert \
	--execute summarize_motif_hits.ipynb \
	--to HTML --output $3/${dtype}_tfm_results \
	--ExecutePreprocessor.timeout=-1 &

wait

python break_overlaps_with_cwm_score.py -o $TFM_MOODS_DIR  -tfm $TFM_TFM_PATH -bw $TFM_SHAP_PATH -g $TFM_REFERENCE_PATH

