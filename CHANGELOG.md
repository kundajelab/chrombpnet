# Changelog

## 2.0.0.dev0 (unreleased)

ChromBPNet moves from TensorFlow 2.8 / tf.keras to Keras 3 on the JAX backend, so it runs natively on CUDA 13 GPUs
(H100, B200, RTX PRO 6000 Blackwell). Commands, inputs and output files are unchanged except where noted below.

### Platform and packaging
- Keras 3 (>= 3.15.1) on JAX (>= 0.11.2). TensorFlow is no longer installed or imported. GPU support comes from
  `jax[cuda13]` (NVIDIA driver >= 580) or the `cuda12` fallback. The CUDA libraries are pip wheels: no system
  CUDA, no `module load cuda`, no `LD_LIBRARY_PATH`.
- Python >= 3.12 (3.13 in the lock files), NumPy 2, pandas 3, pyfaidx >= 0.8.1, weasyprint >= 66.
- `pyproject.toml` replaces `setup.py` and `requirements.txt`. pixi (`pixi.lock`) provides the environments
  `default`, `dev`, `cuda13`, `cuda13-dev` and `cuda12`, including bedtools, bedGraphToBigWig, samtools, MEME and
  pango. uv (`uv.lock`) provides the extras `cuda13` and `cuda12`.
- Removed dependencies: tensorflow, tensorflow-probability, protobuf, deepdish, deeplift, kundajelab-shap,
  modisco-lite 2.0.7, modisco 0.5.16 and matplotlib-venn.
- `modisco` >= 2.5.2 (tfmodisco-lite) replaces modisco-lite. The reports call `modisco report-simple`, because
  `modisco report` writes a different report in 2.5.
- The Docker image is rebuilt on debian:trixie-slim with pixi (`Dockerfile`; CUDA from the JAX wheels). The
  Google Cloud SDK (`gsutil`/`gcloud`) that the 1.x image installed is now opt-in:
  `docker build --build-arg WITH_GCLOUD=true`. The image ignores a host `PYTHONPATH`/`PYTHONHOME` and the
  `~/.local` user site (`PYTHONNOUSERSITE=1`), which Apptainer would otherwise pass in; the Apptainer examples
  also use `--cleanenv`. New GitHub Actions CI checks the lock files, runs the CPU tests on Linux and macOS, runs
  lint, and tests a uv-only install.
- `.github/workflows/docker.yml` replaces the Docker Hub workflow (`latest-docker-image.yml`, which pushed
  `kundajelab/chrombpnet:latest`). On `v*` tags it pushes `<version>-cuda13`, `latest-cuda13` and
  `<sha>-cuda13` images, and the same with `-cuda12`, to the repository's GitHub Container Registry
  (`ghcr.io/kundajelab/chrombpnet`); it can also be run by hand (`workflow_dispatch`). It needs no secrets.
- A job template for SLURM GPU clusters (`workflows/slurm/`).

### Interpretation
- A native JAX DeepSHAP replaces kundajelab-shap (TF1 graph mode). The semantics are the same: the DeepLIFT
  rescale rule for ReLUs, 20 dinucleotide-shuffled references per sequence, the same profile-head weighting, and
  hypothetical scores projected onto the input. It is batched on the GPU and faster than 1.x. Precision `auto`
  (default): full float32 on CPU, TF32 on GPU as in 1.x (full-float32 convolutions are much slower on some GPUs,
  e.g. Blackwell RTX PRO 6000; request them with `--shap-precision highest`). With identical references the
  scores match 1.x within the tolerances of the parity tests (see [VALIDATION.md](VALIDATION.md), section 6).
- DeepSHAP references are now seeded (`--shap-seed`, default 1234) from the seed and each sequence's content.
  Scores are reproducible and do not depend on the order or batching of regions. The 1.x references were
  unseeded, so old and new scores agree statistically, not exactly.
- Contribution-score files keep the 1.x layout: `/raw/seq` int8, `/shap/seq` and `/projected_shap/seq` float16,
  all (N, 4, L), Blosc-compressed. h5py + hdf5plugin now write them instead of deepdish.
- DeepSHAP picks its batch size (unless `--shap-batch-seqs` is given) from the model size and the free GPU
  memory. JAX releases nothing until the process ends, so in `chrombpnet pipeline` that memory stays reserved
  through the TF-MoDISco step. On a shared GPU, set `XLA_PYTHON_CLIENT_MEM_FRACTION`.

### Behaviour changes
- EarlyStopping restores the weights of the best epoch even when training runs all `--epochs` without stopping
  early. 1.x kept the last epoch in that case.
- Initial weights and random streams differ from TensorFlow for the same `--seed`, so retrained models are
  statistically, not bitwise, equivalent to 1.x models. The seed now also controls the initialisation of the
  no-bias ChromBPNet model.
- New model files (`bias.h5`, `bias_model_scaled.h5`, `chrombpnet.h5`, `chrombpnet_nobias.h5`) keep their names
  but are written by Keras 3, so chrombpnet 1.x and TensorFlow 2.x cannot load them (`chrombpnet export` writes a
  copy they can load, see below). chrombpnet 1.x models still load in 2.x: on load, a registered `LogSumExp` layer
  replaces the logsumexp `Lambda` of `chrombpnet.h5`.
- JAX allocates GPU memory on demand: importing chrombpnet sets `XLA_PYTHON_CLIENT_PREALLOCATE=false` unless you
  set it yourself. The pixi environments (so `gpu-check` too) and the Docker image also set it; under `pixi run`
  the environment's value wins over an `export`.
- `bpnet_model.py` and `chrombpnet_with_bias_model.py` are no longer installed on `PATH`. `-a` defaults to the
  packaged files, so only scripts that pass `-a $(which chrombpnet_with_bias_model.py)` need changing. Locate the
  files with `python -c 'import chrombpnet.training.models.chrombpnet_with_bias_model as m; print(m.__file__)'`
  (or `...models.bpnet_model`).
- `reformat_chrombpnet_h5` writes only `chrombpnet_recompiled.h5`. The TensorFlow SavedModel export
  (`chrombpnet_recompiled/`) was dropped.
- `find_chrombpnet_hyperparams` draws the nonpeak subsample for the outlier thresholds with `--seed` (in
  `chrombpnet pipeline` / `train`; run standalone, it has no `--seed` and stays unseeded). The thresholds
  (`counts_sum_min/max_thresh`), `counts_loss_weight` and the filtered peak/nonpeak BEDs are reproducible, and
  they differ from those of any single 1.x run.
- Counts-head DeepSHAP of a full `chrombpnet.h5` is refused, because its counts output is a logsumexp of the
  bias and no-bias counts. 1.x wrote `counts_scores.h5` for it. Use `chrombpnet_nobias.h5` (as the pipelines do),
  or `-pc profile`.
- The reads used for Tn5/DNase shift estimation are a seeded uniform sample instead of an unseeded `shuf -n`.
  The pipelines use the fixed seed 1234, and `-s/--seed` sets it on `reads_to_bigwig` and `auto_shift_detect`.
  GNU `shuf` is no longer needed.
- Gzipped fragment/tagAlign inputs are decompressed in Python instead of by `zcat`. A truncated or corrupt `.gz`
  now fails loudly. It used to silently yield a shift estimate and bigWig built from the reads before the damage.
- TF-MoDISco for the profile and counts heads of `bias pipeline` / `bias qc` still runs one head after the
  other, as in 1.x. Its numba threads follow `NUMBA_NUM_THREADS` if set, otherwise the job's CPUs
  (`SLURM_CPUS_PER_TASK`, else the CPU affinity mask).
- Removed: the legacy modisco-0.5 scripts (`evaluation/modisco/{run_modisco,fetch_tomtom,visualize_motif_matches}.py`,
  `modisco.sh`) and `evaluation/invivo_footprints/`.

### New command
- `chrombpnet export -m model.(h5|keras) -o out.h5 [--legacy-h5] [--count-head {bytecode,named}]` writes a bias,
  chrombpnet or chrombpnet_nobias model as a TF-Keras 2.x full-model .h5 file, in the layout chrombpnet 1.x wrote
  (TF-Keras 2.12): Keras 2 `model_config`, `model_weights/<layer>/<layer>/kernel:0` datasets with the
  `layer_names` / `weight_names` attributes, nested bias / no-bias models stored as in 1.x `chrombpnet.h5`. Keras 3
  auto-generated names are written as TF-Keras named them (`functional` -> `model`, `add_4` ... -> `add` ...).
  TF-Keras 2.x (tested with TF 2.8 and 2.12, `load_model(..., compile=False)`) and bpnet-lite's
  `BPNet.from_chrombpnet` / `ChromBPNet.from_chrombpnet` read these files; 1.x files re-exported this way have the
  same datasets, attributes and model config (byte for byte). There is no `training_config` or optimizer state. The count head of a full chrombpnet model is
  the logsumexp `Lambda` of the 1.x files (Python 3.8 bytecode), so TF-Keras on Python 3.8 - 3.10, and the
  `load_model_wrapper` of chrombpnet 1.x and the variant-scorer, load it without extra custom objects.
  `--count-head named` is for TF-Keras on Python >= 3.11, which cannot unmarshal that bytecode: the `Lambda` then
  names its function and needs `custom_objects={"chrombpnet_logsumexp": ...}` (see
  `chrombpnet/helpers/postprocessing/README.md`). `model_io.load_model` reads both and maps them to `LogSumExp`.

### New optional flags (the defaults keep the 1.x behaviour)
- Training: `--optimizer {adam,muon}`, `--muon-lr`, `--ema`, `--lr-schedule {constant,cosine}`,
  `--precision {default,highest,bf16}`, `--device {auto,gpu,cpu}`.
- Interpretation: `--interpret-subsample` (default 30000), `--shap-seed` (default 1234), `--shap-batch-seqs`,
  `--shap-precision {auto,highest,default}`.
- MoDISco: `--modisco-max-seqlets` (default 50000), `--modisco-window` (default 500), `--tomtom-lite`.
- Preprocessing helpers: `-s/--seed` (default 1234) on `python -m chrombpnet.helpers.preprocessing.reads_to_bigwig`
  and `... auto_shift_detect`, for the reads sampled for shift estimation.

### Validation
- [VALIDATION.md](VALIDATION.md): on ENCODE ENCSR763OVZ (ATAC-seq), ENCODE's 1.x models load and reproduce its
  predicted-signal bigWigs (median per-region r 1.00000), DeepSHAP agrees with ENCODE's scores as closely as two
  of our own runs with different reference seeds, and fold 0 retrained with the default Adam recipe reaches the
  test counts Pearson of ENCODE's fold-0 model (section 5).

### Known issues (pre-existing in 1.x)
- `chrombpnet qc` and `chrombpnet bias qc` read `auxiliary/filtered*.bed` (e.g. `filtered.peaks.bed`,
  `filtered.bias_peaks.bed`) from the output directory they create, so they fail with `FileNotFoundError`
  unless those files are copied there from the training run.

##  Version - 1.5
- Fixed issue #150, regions_used not found while generating bigwigs from impotance h5s

## Version - 1.4
- (MAJOR) Bug in chrombpnet modisco_motifs command. seqlets was limited to 50000. If users wanted to change it to 1 million this did not happen.
- Filter peaks at edges for pred_bw command and bias pipleline. So bias evaluation now done on these filtered peaks.
- Preprocessing deafulted to use unix sort. Provided option to switch to bedtools sort.
- Provided option to use filter chromosomes option in preprocessing.

## Version - 1.3 - Inworks - 2022-12-11
- Changed pipelines to use modisco-lite, old modisco will soon be removed
- Added automatic shift scripts to repo and integrated with the pipeline
- modisco now outputs both html and pdf. pdf can be shared with anyone.
- Use ATAC, DNASE not ATAC_PE, DNASE_SE anymore.
- Simplyfing workflows to include only two main workflows, chrombpnet_train_tf_model, chrombpnet_train_bias_model
- Restructuring README, moving tutorial to additional documentation and introducing FAQ, and only two pipelines (chrombpnet_train_tf_model, chrombpnet_train_bias_model)

## [Unreleased] - 2022-02-28
- Typo fix - (PR#31-36)
- We can now specify ylimit in marginal footprinting plots (PR#27)
- PR#30 merged  to do fast genome wide gc binning and bug fix to ensure case-sensitve GC calculation. Adds a unit test to make checks.

## [Unreleased] - 2022-02-11
- Updated modisco version - (PR#25)

## [Unreleased] - 2022-02-07
- Typo fix - (PR#23-24)

## [Unreleased] - 2022-01-31
- Marginal footprinting no longer hardcoded to take only chr1, it now inputs the fold json. (PR#20)
- Invivo footprinting is made compatible with the new input repo changes. (PR#22)
- Pandas read bugs fixed for gc-matching scripts. (PR#15)
- Tutorial now exits when one part of the pipeline breaks instead of proceeding further. (PR#18)
- Added the fast binning scripts at src/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets (PR#11)
- Added pseudocount in metrics.py (PR#13) - needs further review.

## [Unreleased] - 2022-01-24
- Seed setting fixed for training. Seed was set improperly - preventing shuffling for every epoch. (PR#19)
- return coordinates while training for debugging - returns status of the point such as revcomp,peak/nonpeak,coordinates etc

## [Unreleased] - 2022-01-15
- In `step6_train_chrombpnet_model` changed the file saving name bug  
- Changed scripts in `src/evaluatuion/make_bigwigs/` to be compatible with new dir structure

## [Unreleased] - 2022-01-13
- Added a note in README that users can continue to use the pre-trained bias model from chrombpnet-lite repo, if thats what they have been using till now
- Changed the default modisco crop setting in tutorial from 1000 to 500
- Fixed a bug in `src/evaluation/modisco/run_modisco.py` in directory creation
