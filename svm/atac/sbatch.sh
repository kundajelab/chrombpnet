for task in GM12878 HEPG2 H1ESC IMR90
do
    for fold in `seq 0 9` 
    do
	sbatch -J $task.$fold -o logs/$task.$fold.o -e logs/$task.$fold.e --mem=30000 --time=1440 --mincpus=16 -p akundaje,euan,normal train_predict.sh $task $fold 
    done
done
