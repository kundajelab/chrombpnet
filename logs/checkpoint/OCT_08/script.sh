gpu=1
neg_bed_test=$PWD/results/chrombpnet/ATAC/K562/negatives_data/bpnet.inputs.test.0.negatives.bed
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
model_dir=$PWD/results/chrombpnet/ATAC/K562/4_4_shifted_ATAC_09.12.2021_bias_filters_500/final_model_step3/unplug/
vset=set_1

CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 10 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 25 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 50 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 100 --model_dir $model_dir --motif_type $vset

vset=set_2
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 10 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 25 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 50 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 100 --model_dir $model_dir --motif_type $vset
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 5 --model_dir $model_dir --motif_type $vset

vset=set_1
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/combinatorial_footprints.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --dist 5 --model_dir $model_dir --motif_type $vset






