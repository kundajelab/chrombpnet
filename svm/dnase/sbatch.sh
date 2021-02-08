for task in GM12878 H1ESC HEPG2 K562 IMR90
do
    for fold in `seq 0 9` 
    do
	sbatch -J $task.$fold -o logs/$task.$fold.o -e logs/$task.$fold.e --mem=30000 --time=1440 --mincpus=16 -p akundaje,euan,owners train_predict.sh $task $fold 
	#sbatch -J $task.$fold -o logs/$task.$fold.o -e logs/$task.$fold.e --mem=30000 --time=1440 --mincpus=16 -p akundaje,euan,owners train_predict_gkm_rbf.sh $task $fold 
    done
done
