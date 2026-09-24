
# Scripts to get contribution scores from ChromBPNet models

The scripts in this folder are used to do get sequence contribution scores for the chrombpnet models. They run a native JAX implementation of DeepSHAP (`explainer.py`, `deeplift_jax.py`) that reproduces the chrombpnet 1.x `shap` (kundajelab TFDeepExplainer) scores; `shap`, `deeplift` and TensorFlow are no longer needed. Chrombpnet models have two heads - (1) profile head (to predict the shape of a profile) and (2) counts head  (to predict the total counts in a profile). The scripts can provide sequence contribution scores for both the heads.

## Usage

```
python -m chrombpnet.evaluation.interpret.interpret -g [genome_fasta] -r [regions] -m [model_h5] -o [output_prefix] -d [debug_chr] -p [profile_or_counts]
```

`chrombpnet contribs_bw` runs the same step and then writes the bigwigs.

## Example Usage

```
python -m chrombpnet.evaluation.interpret.interpret -g /path/to/genome.fasta -r /path/to/regions.bed -m /path/to/model.h5 -o /path/to/output_dir/output_prefix -p counts profile
```

## Input Format
- genome_fasta: Reference geneome fasta.
- regions: Must be a narrowPeak file, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region is centered at start + summit internally, across all regions.
- model_h5: Model in hdf5 format.
- output_prefix: Output prefix path and name to use. See Output format section below to understand how this is used.
- profile_or_counts: Provide either only `profile` or `counts` or both with space speration. The default setting for this uses both `profile` and `counts`. This argument is used to indicate the output head (counts/profiles or both) to use for getting contribution scores. 
- `--seed` (default 1234): seed of the 20 dinucleotide-shuffled references per region. Each region's references are seeded from the seed and the region's sequence, so scores are reproducible and do not depend on region order or batching. chrombpnet 1.x references were unseeded, so old and new scores agree only up to shuffle noise.
- `--batch-seqs`: regions per device step (default: from the model's activation size and the free device memory, with steps of at most ~12 GB); halved automatically when the GPU runs out of memory.
- `--precision` (default `auto`): full float32 on CPU, TF32 on GPU (as chrombpnet 1.x did). `highest` forces full float32 convolutions on the GPU too, which is much slower on some GPUs (e.g. an RTX PRO 6000 Blackwell); `default` always uses the backend default. With identical references, full-float32 scores match 1.x within the tolerances of `tests/parity/test_shap_parity.py` (see [VALIDATION.md](../../../VALIDATION.md), section 6), and TF32 scores agree closely with the float32 ones.
- `--device` (`auto`, `gpu`, `cpu`), `--save-references` (writes `output_prefix`.references.npz), `--references` (reuse such a file), `--profile-weighting` (`chrombpnet` reproduces 1.x; `softmax_x`, `tangermeme` and `chrombpnet_tf` are for comparisons only).
- The counts head of a full `chrombpnet.h5` (logsumexp of the bias and no-bias counts) is not supported: interpret `chrombpnet_nobias.h5`, as the pipeline does. Its profile head is supported.
- debug_chr: This argument is optional and can be used for debugging if needed. Mention a chromsome name to get contribution scores only for regions in that chromosome. The chromsome name should be present in the first column of the `regions` bed.

## Output Format

The following two files are created using the `output_prefix` as prefix for the output when both `profile_or_counts` includes both `profile` and `counts`. If only one of the option is used, only one of the files is generated. Note that prefix can include a directory path and name for the output file. 


- `output_prefix`.profile_scores.h5: A h5 file containing a dictionary with contribution scores for the profile head of the chrombpnet model. The dictionary consists of the following keys `raw`, `shap` and `projected_shap`. The values of each of these keys are explained below. All the values have the shape `(L,4,inputlen)` where  L is the number of inputs given in `regions` bed file and `inputlen` is inferred from the input length of the `model_h5`.
    - `raw` key consists of  one hot encoding of a sequence from the `regions` bed file.
    - `shap` key consists of hypothetical contribution scores obtained from doing deepshap. These values will be used in the MODISCO pipeline (for de-novo motif discovery) next.
    - `projected_shap` key consists of values obtained by multiplying `shap` values with `raw` values. These values will be used in building browser tracks.

- `output_prefix`.counts_scores.h5: A h5 file containing a dictionary with contribution scores for the count head of the chrombpnet model. Has similar structure as the `output_prefix`.profile_scores.h5.

- The h5 files use the Blosc filter (as chrombpnet 1.x did through deepdish): `import hdf5plugin` before opening them with h5py. `/raw/seq` is int8, `/shap/seq` and `/projected_shap/seq` are float16.

- `output_prefix`.interpret.args.json: the arguments and the DeepSHAP settings used (seed, precision, device, batch size, JAX devices).

- `output_prefix`.interpreted_regions.bed: A bed file with regions used for interpretations. The script filters regions whose resulting sequences dont have the same length as `input_len` inferred from the model.  This bed file might have regions fewer than `regions` provided as input. 


