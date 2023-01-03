# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
reference_fasta=$1
regions=$2
model_h5=$3
output_dir=$4

## deepshap run

python $PWD/src/evaluation/interpret/interpret.py \
       --genome=$reference_fasta \
       --regions=$regions \
       --output_prefix=$output_dir/bias \
       --model_h5=$model_h5 \

## modisco run

python src/evaluation/modisco/run_modisco.py -s $output_dir/bias -p "profile" -o $output_dir -m 50000 -c 500
python src/evaluation/modisco/run_modisco.py -s $output_dir/bias -p "counts" -o $output_dir -m 50000 -c 500
