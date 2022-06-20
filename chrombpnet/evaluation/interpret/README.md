
# Scripts to get contribution scores from ChromBPNet models

The scripts in this folder are used to do get sequence contribution scores for the chrombpnet models. These scripts require the `shap` installation. Chrombpnet models have two heads - (1) profile head (to predict the shape of a profile) and (2) counts head  (to predict the total counts in a profile). The scripts can provide sequence contribution scores for both the heads.

## Usage

```
chrombpne_deepshap -g [genome_fasta] -r [regions] -m [model_h5] -o [output_prefix] -d [debug_chr] -p [profile_or_counts]
```

## Example Usage

```
chrombpnet_deepshap -g /path/to/genome.fasta -r /path/to/regions.bed -m /path/to/model.h5 -o /path/to/output_dir/output_prefix -p counts profiles
```

## Input Format
- genome_fasta: Reference geneome fasta.
- regions: Must be a narrowPeak file, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region is centered at start + summit internally, across all regions.
- model_h5: Model in hdf5 format.
- output_prefix: Output prefix path and name to use. See Output format section below to understand how this is used.
- profile_or_counts: Provide either only `profile` or `counts` or both with space speration. The default setting for this uses both `profile` and `counts`. This argument is used to indicate the output head (counts/profiles or both) to use for getting contribution scores. 
- debug_chr: This argument is optional and can be used for debugging if needed. Mention a chromsome name to get contribution scores only for regions in that chromosome. The chromsome name should be present in the first column of the `regions` bed.

## Output Format

The following two files are created using the `output_prefix` as prefix for the output when both `profile_or_counts` includes both `profile` and `counts`. If only one of the option is used, only one of the files is generated. Note that prefix can include a directory path and name for the output file. 


- `output_prefix`.profile_scores.h5: A h5 file containing a dictionary with contribution scores for the profile head of the chrombpnet model. The dictionary consists of the following keys `raw`, `shap` and `projected_shap`. The values of each of these keys are explained below. All the values have the shape `(L,4,inputlen)` where  L is the number of inputs given in `regions` bed file and `inputlen` is inferred from the input length of the `model_h5`.
    - `raw` key consists of  one hot encoding of a sequence from the `regions` bed file.
    - `shap` key consists of hypothetical contribution scores obtained from doing deepshap. These values will be used in the MODISCO pipeline (for de-novo motif discovery) next.
    - `projected_shap` key consists of values obtained by multiplying `shap` values with `raw` values. These values will be used in building browser tracks.

- `output_prefix`.count_scores.h5: A h5 file containing a dictionary with contribution scores for the count head of the chrombpnet model. Has similar structure as the `output_prefix`.profile_scores.h5.

- `output_prefix`.interpreted_regions.bed: A bed file with regions used for interpretations. The script filters regions whose resulting sequences dont have the same length as `input_len` inferred from the model.  This bed file might have regions fewer than `regions` provided as input. 


