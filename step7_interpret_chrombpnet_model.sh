reference_fasta=$1
regions=$2
model_h5=$3
output_dir=$4
codedir=$5

## deepshap run

python $codedir/src/evaluation/interpret/interpret.py \
       --genome=$reference_fasta \
       --regions=$regions \
       --output_prefix=$output_dir/corrected \
       --model_h5=$model_h5 \

## modisco run

python $codedir/src/evaluation/modisco/run_modisco.py -s $output_dir/corrected -p "profile" -o $output_dir -m 50000 -c 1000
python $codedir/src/evaluation/modisco/run_modisco.py -s $output_dir/corrected -p "counts" -o $output_dir -m 50000 -c 1000
