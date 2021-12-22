reference_fasta=$1
regions=$2
model_h5=$3
output_dir=$4


    
#python $PWD/src/evaluation/interpret/interpret.py \
#        --genome=$reference_fasta \
#        --regions=$regions \
#        --output_prefix=$output_dir/interpret \
#        --model_h5=$model_h5 \

## modisco run

python src/evaluation/modisco/run_modisco.py -s $output_dir/interpret -p "profile" -o $output_dir -m 50000 -c 1000
python src/evaluation/modisco/run_modisco.py -s $output_dir/interpret -p "counts" -o $output_dir -m 50000 -c 1000
