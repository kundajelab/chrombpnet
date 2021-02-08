prefix=/oak/stanford/groups/akundaje/projects/enzymatic_bias_correction/svm
for cell_line in K562 HEPG2
do
    for split in `seq 0 4`
    do
	for chunk in `seq 0 9`
	do
	    cat $prefix/pred/gkmexplain.$cell_line.$split.$chunk.dnase.txt >> $prefix/aggregate/gkmexplain.$cell_line.$split.txt 
	    
	done
    done
done
