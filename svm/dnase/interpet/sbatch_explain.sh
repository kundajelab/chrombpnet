for cell_line in GM12878  H1ESC IMR90
do
    for split in `seq 0 9`
    do
	sbatch -J $cell_line.$split -o logs_explain/$cell_line.$split.o -e logs_explain/$cell_line.$split.e -t 1-0 -p akundaje,euan,owners explain.sh 
    done
done
