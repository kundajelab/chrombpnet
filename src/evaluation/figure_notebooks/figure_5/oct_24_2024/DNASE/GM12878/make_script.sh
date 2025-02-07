

celltype="GM12878_new"
mode="counts"
dir='/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/'$celltype'/merge_folds_new_may_05_24/'$mode
modisco convert-backward -i $dir/modisco_$mode".h5" -o $dir/modisco_old_format.h5



mode="profile"
dir='/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/'$celltype'/merge_folds_new_may_05_24/'$mode
modisco convert-backward -i $dir/modisco_$mode".h5" -o $dir/modisco_old_format.h5



