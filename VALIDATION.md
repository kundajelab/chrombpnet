# Validation of ChromBPNet 2.x (Keras 3 / JAX)

This page compares ChromBPNet 2.x with chrombpnet 1.x (TensorFlow 2.x) on one public ENCODE ATAC-seq dataset,
using the models and bigWigs that ENCODE released, and describes how the code is tested. Sections 2 to 5 are the
comparisons, section 6 is the test suite, sections 7 and 8 are GPU notes and what has not been tested yet.

## 1. What was compared, and on what

The 1.x side of every comparison in sections 2 to 5 is the files ENCODE released. chrombpnet 2.x does not
install TensorFlow, and none was installed for these runs.

| | |
|---|---|
| Experiment | [ENCSR763OVZ](https://www.encodeproject.org/ENCSR763OVZ/) (ATAC-seq) |
| ChromBPNet model annotation | [ENCSR729FNL](https://www.encodeproject.org/ENCSR729FNL/) |

| ENCODE file | content | used in section |
|---|---|---|
| [ENCFF352NED](https://www.encodeproject.org/files/ENCFF352NED/) | models tar: 5 folds, each with `chrombpnet`, `chrombpnet_nobias` and `bias_scaled` (1.x `.h5`) | 2, 3, 4; 5 (fold-0 bias model) |
| [ENCFF867HHL](https://www.encodeproject.org/files/ENCFF867HHL/) | training / test regions tar | 5 |
| [ENCFF683KQI](https://www.encodeproject.org/files/ENCFF683KQI/) | regions selected for the predicted-signal and contribution-score bigWigs | 3, 4 |
| [ENCFF477SDP](https://www.encodeproject.org/files/ENCFF477SDP/) | predicted signal bigWig, full model (`chrombpnet`), 5-fold mean | 3 |
| [ENCFF173FZC](https://www.encodeproject.org/files/ENCFF173FZC/) | bias-corrected predicted signal bigWig (`chrombpnet_nobias`), 5-fold mean | 3 |
| [ENCFF783VNM](https://www.encodeproject.org/files/ENCFF783VNM/) | counts contribution-score bigWig, 5-fold mean | 4 |
| [ENCFF955CFK](https://www.encodeproject.org/files/ENCFF955CFK/) | profile contribution-score bigWig, 5-fold mean | 4 |
| [ENCFF185BWD](https://www.encodeproject.org/files/ENCFF185BWD/) | observed signal bigWig that ENCODE's models were trained on | 5 |

All runs used one machine:

| | |
|---|---|
| GPU | NVIDIA RTX PRO 6000 Blackwell Server Edition (sm_120, 96 GB) |
| Driver | 595.71.05 (CUDA 13.2) |
| Software | chrombpnet 2.0.0.dev0, keras 3.15.1, jax 0.11.2 (`jax[cuda13]`) |
| Precision | default: TF32 convolutions on the GPU, as TensorFlow used for 1.x |

## 2. Loading 1.x models

All 15 released models (5 folds x 3 model types) load with `chrombpnet.training.utils.model_io.load_model`,
without TensorFlow. The loader replaces the logsumexp `Lambda` count head of `chrombpnet.h5` with a registered
`LogSumExp` layer and does not run the Python bytecode stored in the file.

| model | `keras_version` in the file | structure | files loaded |
|---|---|---|---|
| `bias_scaled` | 2.4.0 | bias model | 5 / 5 |
| `chrombpnet` | 2.8.0 | nested legacy models (bias + no-bias), `Lambda` logsumexp count head | 5 / 5 |
| `chrombpnet_nobias` | 2.4.0 | no-bias model | 5 / 5 |

As a consistency check, the full model's outputs were compared with its two parts computed separately: profile
logits of `chrombpnet` = logits of `chrombpnet_nobias` + logits of `bias_scaled`, and logcounts of `chrombpnet` =
logsumexp of the two logcounts.

| output | max absolute difference, full model vs parts |
|---|---|
| logits | 1.91e-06 |
| logcounts | <= 4.77e-07 |

## 3. Predictions vs ENCODE

- **Regions:** 300 isolated summits on the fold-0 test chromosomes (chr1, chr3, chr6), sampled with seed 1234.
- **Predictions:** each of the 5 released folds, averaged the way ENCODE did for its bigWigs:
  `softmax(mean logits) * (exp(mean logcounts) - 1)`.
- **Comparison:** per region, over the 1000 bp output window, against the values in ENCODE's predicted-signal
  bigWig for that model.

| model | ENCODE bigWig | median Pearson r | 5th percentile r | median ratio of summed signal (5th-95th pct.) | median max. abs. error / region max. |
|---|---|---|---|---|---|
| `chrombpnet` | ENCFF477SDP | 1.00000 | 1.00000 | 1.0000 (0.9999-1.0000) | 0.0003 |
| `chrombpnet_nobias` | ENCFF173FZC | 1.00000 | 1.00000 | 0.9999 | 0.0003 |

The predicted tracks match ENCODE's to five decimal places of correlation. In the median region, the largest
difference is 0.03% of the region's maximum.

## 4. Contribution scores vs ENCODE

- **Method:** the native JAX DeepSHAP of 2.x on the `chrombpnet_nobias` model, counts and profile heads, 20
  dinucleotide-shuffled references per region, 2114 bp inputs, run for each of the 5 folds and averaged.
- **Comparison:** the 5-fold mean against ENCODE's 5-fold mean bigWigs (ENCFF783VNM counts, ENCFF955CFK
  profile), on the same 300 regions as section 3, per-region Pearson r over the central summit +/- 500 bp.

| head | vs ENCODE: median r | vs ENCODE: 5th percentile r | our run vs our run, two shuffle seeds: median r |
|---|---|---|---|
| counts | 0.9699 | 0.9600 | 0.9696 |
| profile | 0.9318 | 0.8365 | 0.9316 |

**Why the correlation is not 1: the noise floor.** DeepSHAP scores depend on the randomly shuffled reference
sequences. ENCODE's 1.x scores used unseeded dinucleotide shuffles, so their references cannot be regenerated.
To measure how much of the difference the references alone explain, DeepSHAP was run twice on the same models
with two different shuffle seeds, and the two runs were compared in the same way (last column). Our agreement
with ENCODE (0.9699 counts, 0.9318 profile) equals our agreement with ourselves (0.9696, 0.9316). The difference
from ENCODE is the size expected from the choice of references, not a difference introduced by the port. The
profile head is more sensitive to the references than the counts head in both comparisons.

- **Additivity (completeness):** for every region, the scores sum to the change in the explained output between
  the sequence and its references, with a maximum relative error <= 0.0002.
- **Throughput:** about 10 regions/s per fold (counts and profile heads, 20 references per region, 2114 bp).

## 5. Retraining from ENCODE's bigWig: Adam vs Muon

Fold 0 was retrained from ENCODE's observed-signal bigWig (ENCFF185BWD), with ENCODE's fold-0 regions and
ENCODE's fold-0 bias model, using the 1.x recipe:

| setting | value |
|---|---|
| filters / dilated layers | 512 / 8 |
| input / output length | 2114 / 1000 bp |
| batch size, learning rate | 64, 1e-3 |
| epochs | at most 50, early stopping patience 5 (the best epoch is restored) |
| negative sampling ratio | 0.1 |
| filtered peaks, train / valid / test | 27,591 / 2,821 / 8,851 |
| `counts_loss_weight` | 3.7 (the same as ENCODE's model) |

Adam at 1e-3 is the 1.x recipe and the default. Muon (`--optimizer muon`, default `--muon-lr 2e-3`) is opt-in:
it updates the dilated convolution kernels, and Adam updates the other weights. Everything else was identical.
Each optimizer was run with three seeds. Metrics are on the test peaks (fold-0 test chromosomes chr1, chr3,
chr6).

| optimizer | seed | best val_loss @ epoch / epochs run | test counts Pearson | test counts Spearman | median JSD | s/epoch |
|---|---|---|---|---|---|---|
| adam | 1234 | 136.745 @ 7 / 12 | 0.5996 | 0.6833 | 0.6947 | n/a (timing not comparable: GPU shared) |
| adam | 1 | 136.516 @ 6 / 11 | 0.5962 | 0.6833 | 0.6938 | 46.4 |
| adam | 2 | 136.589 @ 5 / 10 | 0.5992 | 0.6831 | 0.6926 | 46.4 |
| muon (lr 2e-3) | 1234 | 136.586 @ 3 / 8 | 0.6038 | 0.6856 | 0.6937 | 47.2 |
| muon (lr 2e-3) | 1 | 136.800 @ 3 / 8 | 0.6041 | 0.6838 | 0.6947 | 47.2 |
| muon (lr 2e-3) | 2 | 136.807 @ 2 / 7 | 0.5945 | 0.6759 | 0.6944 | 47.2 |

Means over the 3 seeds (+/- standard deviation):

| optimizer | best val_loss | best epoch | steps to best | test counts Pearson | test counts Spearman | median JSD |
|---|---|---|---|---|---|---|
| adam | 136.616 +/- 0.117 | 6.0 +/- 1.0 | 2375-3325 | 0.5983 +/- 0.0018 | 0.6832 | 0.6937 |
| muon | 136.731 +/- 0.126 | 2.7 +/- 0.6 | 950-1425 | 0.6008 +/- 0.0054 | 0.6818 | 0.6943 |

The cost per step is the same for both (about 95 ms per step at batch 64).

**Compared with ENCODE's released fold-0 model (1.x):** test counts Pearson 0.598, Spearman 0.684. The Adam
retraining gives 0.5983 and 0.6832. JSD is not compared with ENCODE, because ENCODE's QC report and chrombpnet's
metrics define it differently.

**Conclusion.** With the default Adam recipe, 2.x retrains a model with the same test counts correlation as
ENCODE's 1.x model. Muon reached its best validation loss in about half the epochs (2.7 vs 6.0, 950-1425 vs
2375-3325 steps), with the same final accuracy within seed noise: its mean best validation loss is slightly
higher (136.731 vs 136.616) and its mean test Pearson slightly higher (0.6008 vs 0.5983), differences about the
size of the seed-to-seed standard deviations (0.117-0.126 and 0.0018-0.0054). This is 3 seeds on one dataset
and one fold. Adam stays the default, and Muon stays an opt-in.

## 6. How the code is tested

### Unit tests (CPU, in CI)

The 19 test modules at the top level of `tests/` (about 210 test functions) run on the CPU. `tests/conftest.py` builds a small
synthetic dataset once per session: 3 chromosomes (~90 kb) with a FASTA and `.fai`, a chrom sizes file,
10-column narrowPeak peaks and nonpeaks, a counts bigWig with read pile-ups at the peak summits, an AP-1 motif
planted at every summit, and a fold JSON. The tests need no network access and no reference data.

| CI job (`.github/workflows/ci.yml`) | what it runs |
|---|---|
| lock files up to date | `pixi lock --check`, `uv lock --check` |
| lint | `pixi run -e dev lint` (ruff) |
| test (ubuntu-latest, macos-14) | `pixi run -e dev test`, i.e. `pytest -q -m 'not gpu and not slow'`, with `JAX_PLATFORMS=cpu` |
| uv only (no conda tools) | `uv sync --locked`, a check that importing chrombpnet selects the JAX backend, and `pytest -q -m "not gpu and not slow and not needs_cli"` |
| cuda13 environment installs (pushes only, not pull requests) | installs the `cuda13` environment on a CPU runner and imports it |

Some of what the unit tests check:

| module | checks |
|---|---|
| `test_losses.py` | `multinomial_nll` against a float64 reference of TensorFlow Probability's `Multinomial(...).log_prob` |
| `test_dinuc_shuffle.py` | the dinucleotide shuffle is bit-identical to deeplift 0.6.13 (a verbatim copy is in the test); references are seeded by sequence content |
| `test_deeplift_jax.py` | DeepLIFT rescale rule, the midpoint rule for the profile head, summation to delta, agreement with a NumPy brute-force DeepLIFT, batching and padding invariance, counts head of a full model refused |
| `test_model_io.py` | loading hand-built TF-Keras 2.x `.h5` files (nested models, the logsumexp `Lambda`), re-saving them with Keras 3, other `Lambda` layers refused |
| `test_export_legacy_h5.py` | `chrombpnet export --legacy-h5`: file layout and Keras 2 model config, round trip through `model_io.load_model` (atol 1e-6), bpnet-lite's reader; with `CHROMBPNET_GOLDENS` set, the layout against the TF-Keras 2.12 files and re-exported 1.x files identical to the originals; with `CHROMBPNET_TF_PYTHON` (a Python with TensorFlow 2.x) set, TF-Keras loads the exports |
| `test_models.py`, `test_training.py`, `test_hyperparams.py`, `test_predict.py`, `test_interpret.py`, `test_footprints.py` | architectures, `train`, hyperparameter search, prediction and bigWig writing, interpretation and footprinting end to end on the synthetic data with tiny models, including the restored best epoch, EMA, Muon, cosine schedule and bf16 options |
| `test_preprocessing.py`, `test_modisco_run.py`, `test_reports.py`, `test_cli.py` | seeded read sampling and shift detection, gzipped inputs, bigWig building, the modisco calls, the HTML/PDF reports, CLI defaults and pipeline wiring |
| `test_imports.py`, `test_infra.py` | every module imports without TensorFlow; pixi activation, CI action references and the Docker image environment |

On a GPU machine, `pixi run -e cuda13-dev test` runs the same tests with the GPU build of JAX, and
`pixi run -e cuda13-dev test-gpu` runs the tests marked `gpu`. At the moment that is one test, which checks that
JAX sees a GPU.

### Parity with 1.x (`tests/parity/`, run by hand)

The goldens are made by `tests/goldens/make_goldens.py` (driven by `tests/goldens/run_goldens.sh`) inside a
container with chrombpnet 1.x (Python 3.8, TensorFlow 2.12, kundajelab-shap), on the CPU. The script imports the
container's own chrombpnet, not this checkout. Its subcommands write the 1.x predictions of a set of models
(`predict`), DeepSHAP with injected, content-seeded references (`shap`), the model outputs on those sequences and
references (`fvals`), and 20 Adam steps from exported initial weights on fixed batches (`trace`).

The parity tests read these environment variables. Each test is skipped when a variable it needs is unset, so
they do not run in CI:

- `CHROMBPNET_GOLDENS`: the directory the goldens were written to.
- `CHROMBPNET_LEGACY_RESULTS`: the chrombpnet 1.x results tree with the models the goldens were made from, laid
  out as `bias_models/bias_model_<b>/<prefix>/{models,evaluation,auxiliary}/` (one directory per bias model) and
  `full_models/<prefix>/{models,auxiliary}/`. `test_modisco_parity.py` also accepts several directories separated
  by `os.pathsep`.
- `CHROMBPNET_LEGACY_PREFIX`: the run name `<prefix>`, which is both the run directory above and the file prefix
  of the bias-model files (`<prefix>_bias.h5`).
- `CHROMBPNET_LEGACY_BIAS`: the bias model `<b>` whose DeepSHAP goldens were written (`shap_bias_<b>*.npz`,
  `fvals_bias_<b>.npz`).
- Optional: `CHROMBPNET_MODEL_BIAS` and `CHROMBPNET_MODEL_NOBIAS`, explicit model paths for
  `test_shap_parity.py`.

`run_goldens.sh` takes the same three values as `RESULTS`, `PREFIX` and `BIAS` (each falls back to the matching
`CHROMBPNET_LEGACY_*` variable), plus `SANDBOX` (the chrombpnet 1.x container) and `GENOME` (the reference FASTA
the models were trained on).

On a GPU, also export `XLA_FLAGS=--xla_gpu_deterministic_ops=true`.

| test file | what it compares | tolerance |
|---|---|---|
| `test_predict_parity.py` | logits and logcounts of each 1.x bias model under `bias_models/` and of the run's `bias_model_scaled`, `chrombpnet` and `chrombpnet_nobias`, at matmul precision `highest` | max. abs. difference <= 1e-4 |
| | full model vs `chrombpnet_nobias` + `bias_model_scaled` (logits added, logcounts by logsumexp) | max. abs. difference <= 1e-4 |
| `test_shap_parity.py` | references regenerated from the sequences | identical to the goldens |
| | the 1.x profile-head multipliers follow the midpoint rule (a check of which rule 1.x applies) | median relative error < 1e-3 |
| | the 1.x counts-head multipliers sum to the output difference | max. relative error < 1e-3 |
| | model outputs on the DeepSHAP inputs, bias model and `chrombpnet_nobias` | atol 1e-4 (CPU) / 1e-3 (GPU) |
| | DeepSHAP multipliers and hypothetical scores, both heads, both models, precision `highest` | relative L2 <= 1e-4 (CPU) / 1e-3 (GPU) |
| `test_training_parity.py` | 20 Adam steps (lr 1e-3, batch 32) from the same initial weights on the same batches, for a 128-filter x 4-layer bias model and a 32 x 4 chrombpnet model with a frozen 1.x bias model: per-step total, profile and counts losses | relative difference < 1e-5 at the first step, < 1e-3 at the last |
| | final weights after the 20 steps | relative L2 < 5e-3 over all weights, < 5e-2 per array (the profile-head output bias is excluded: its true gradient is 0, so its Adam steps follow float noise in both implementations) |
| | a 1.x `bias_model_scaled.h5` used as the pretrained bias stays frozen after a training step | outputs equal at rtol / atol 1e-5 |
| `test_modisco_parity.py` | the report keeps exactly the `motifs.html` rows that the 1.x parser kept, on 1.x reports | same rows |
| | (slow) `modisco motifs` 2.5 rerun on 1.x profile scores vs the 1.x (modisco-lite 2.0.7) patterns | number of patterns within max(3, 30%); at least 80% of the 10 largest 1.x patterns matched by a CWM with Pearson >= 0.9 (offsets up to 5 bp, both strands) |

## 7. GPU notes

While training (512 filters, batch 64) on the RTX PRO 6000, the GPU runs kernels almost all the time, but its
power draw and memory bandwidth stay well below the card's limits: training at batch 64 does not saturate this
GPU. Larger batches and bf16 training are follow-ups. The current `--precision bf16` option loses accuracy and
is not recommended yet. The default precision uses TF32 for convolutions, as TensorFlow did for 1.x. See
[GPU notes](README.md#gpu-notes) in the README for drivers, memory and shared GPUs.

## 8. Limitations and what has not been tested yet

- **Other GPUs.** Only the RTX PRO 6000 Blackwell (sm_120) was used. H100 (sm_90) and B200 (sm_100) have not
  been tested. On such a machine, `pixi run -e cuda13 gpu-check` lists the GPUs JAX sees and
  `pixi run -e cuda13-dev test-gpu` runs the GPU tests.
- **One dataset, one fold.** One ENCODE ATAC-seq dataset was used, and only fold 0 was retrained, with three
  seeds per optimizer and Muon at one learning rate (2e-3). No DNase-seq data was compared.
- **Bias model not retrained.** Retraining reused ENCODE's fold-0 bias model, so bias-model training was not
  compared on this data.
- **Preprocessing not compared.** Retraining started from ENCODE's bigWig, so the steps from reads to bigWig
  (shift detection, bigWig building) are not part of these comparisons.
- **300 regions.** Predictions and contribution scores were compared on 300 regions, and contribution scores
  only for `chrombpnet_nobias`.
- **Downstream outputs.** TF-MoDISco motifs, footprints and the HTML/PDF reports were not compared with ENCODE's.
- **JSD** is not comparable with ENCODE's QC report (different definitions).
- **bf16** loses accuracy and is not recommended; training at batch 64 does not use the whole GPU.
