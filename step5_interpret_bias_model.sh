reference_fasta=$1
regions=$2
model_h5=$3


mkdir output/bias_model/deepshap  
    
#CUDA_VISIBLE_DEVICES=0 python $PWD/src/evaluation/interpret/interpret.py \
#         --genome=$reference_fasta \
#         --regions=$regions \
#         --output_prefix="output/bias_model/deepshap/gm12878" \
#         --model_h5=$model_h5 \


## modisco run

mkdir output/bias_model/deepshap/run_modisco
output_dir="output/bias_model/deepshap/run_modisco"
output_prefix="output/bias_model/deepshap/gm12878"
python src/evaluation/modisco/run_modisco.py -s $output_prefix -p "profile" -o $output_dir -m 50000 -c 1000
python src/evaluation/modisco/run_modisco.py -s $output_prefix -p "counts" -o $output_dir -m 50000 -c 1000
