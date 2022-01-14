reference_fasta=$1
regions=$2
model_h5=$3
output_dir=$4

## deepshap run

python $PWD/evaluation/interpret/interpret.py \
       --genome=$reference_fasta \
       --regions=$regions \
       --output_prefix=$output_dir/bias \
       --model_h5=$model_h5 \

## modisco run

python $PWD/evaluation/modisco/run_modisco.py -s $output_dir/bias -p "profile" -o $output_dir -m 50000 -c 500
python $PWD/evaluation/modisco/run_modisco.py -s $output_dir/bias -p "counts" -o $output_dir -m 50000 -c 500
