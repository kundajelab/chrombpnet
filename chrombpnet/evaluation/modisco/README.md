# TF-MoDISco on ChromBPNet contribution scores

`chrombpnet pipeline` and `chrombpnet bias pipeline` run TF-MoDISco on the contribution scores of a random
subsample of peaks (30,000 by default) and embed the resulting motifs in `evaluation/overall_report.{html,pdf}`:

| pipeline | scores | MoDISco results | report |
|---|---|---|---|
| `chrombpnet pipeline` | `auxiliary/interpret_subsample/chrombpnet_nobias.profile_scores.h5` | `auxiliary/interpret_subsample/modisco_results_profile_scores.h5` | `evaluation/modisco_profile/motifs.html`, `evaluation/chrombpnet_nobias_profile.pdf` |
| `chrombpnet bias pipeline` | `auxiliary/interpret_subsample/bias.{profile,counts}_scores.h5` | `auxiliary/interpret_subsample/modisco_results_{profile,counts}_scores.h5` | `evaluation/modisco_{profile,counts}/motifs.html`, `evaluation/bias_{profile,counts}.pdf` |

MoDISco comes from the PyPI package `modisco>=2.5.2` (tfmodisco-lite; installed with chrombpnet). Do not also
install `modisco-lite` or the legacy `modisco` 0.5: they provide the same module and command.
The two steps are wrapped in `run.py`:

```python
from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.evaluation.modisco.run import modisco_motifs, modisco_report

modisco_motifs("chrombpnet_nobias.profile_scores.h5", "modisco_results_profile_scores.h5",
               max_seqlets=50000, window=500)
modisco_report("modisco_results_profile_scores.h5", "modisco_profile/",
               get_default_data_path(DefaultDataFile.motifs_meme), tomtom_lite=False)
```

which run the equivalent of

```
modisco motifs -i chrombpnet_nobias.profile_scores.h5 -n 50000 -o modisco_results_profile_scores.h5 \
    -w 500 -l 2 -z 20 -f 5 -t 20 -g 5 -j 0
modisco report-simple -i modisco_results_profile_scores.h5 -o modisco_profile/ -m $(print_meme_motif_file) -n 3
```

Both functions also take an optional `threads=` (see [Speed](#speed)).

## Settings

- The pattern settings (`-z 20 -f 5 -t 20 -g 5 -j 0`, 2 Leiden runs) are the ones chrombpnet 1.x used with
  modisco-lite 2.0.7. modisco 2.5 changed the command-line defaults to `-t 30 -g 10` (50 bp patterns), so they
  are always passed explicitly.
- `-n` / `max_seqlets` (`--modisco-max-seqlets`, default 50000) caps the seqlets per metacluster. It is the main
  runtime knob: the all-pairs seqlet similarity step grows about quadratically with it. The cap keeps the first
  seqlets in region order, so when it binds the last regions of the (random) subsample are not used.
- `-w` / `window` (`--modisco-window`, default 500) is the width around the region centre used for motif
  discovery. chrombpnet uses 500 bp of the 2114 bp input to avoid AT-rich nucleosome-flank motifs.
- The report is `modisco report-simple`, which writes `motifs.html`. (`modisco report` in modisco 2.5 is a
  different, descriptive report that writes `report.html`.)

## TOMTOM matches

Each pattern is matched against `chrombpnet/data/motifs.meme.txt` (MEME TF motifs plus recurring Tn5/DNase
bias motifs; `print_meme_motif_file` prints its path) and the top 3 matches are shown.

- Default: MEME's `tomtom` binary (Pearson distance), reporting q-values in columns `qval0..2`, as in
  chrombpnet 1.x. `tomtom` comes from bioconda `meme`, which the linux pixi environments include.
- `--tomtom-lite` (`tomtom_lite=True`): memelite's TOMTOM-lite. No MEME install, and the report step is several
  times faster, but it uses Euclidean distance and reports uncorrected p-values (`pval0..2`), so top matches
  can differ. Use it where MEME is unavailable (e.g. the macOS pixi environment). The HTML reports adapt their
  wording to whichever columns are present.

## Speed

MoDISco runs on the CPU only (numba). `run.py` sets `NUMBA_NUM_THREADS` and `OMP_NUM_THREADS` (to the same
value) for each modisco process, both `modisco motifs` and `modisco report-simple`. The value is chosen in this
order:

1. the `threads=` argument of `modisco_motifs` / `modisco_report`, if given;
2. otherwise `NUMBA_NUM_THREADS`, if you set it to a positive integer (your setting wins over the CPU count);
3. otherwise `OMP_NUM_THREADS`, if you set it to a positive integer;
4. otherwise `SLURM_CPUS_PER_TASK`;
5. otherwise the CPUs this process may run on (the CPU affinity mask, else `os.cpu_count()`).

Steps 2-5 are `run.default_threads()`. To cap MoDISco on a shared machine, export `NUMBA_NUM_THREADS`
before running chrombpnet. `run.py` also points `NUMBA_CACHE_DIR` at a writable directory if it is unset.
The thread count does not change the results. On a GPU node with few CPUs, prefer running the MoDISco step as
its own CPU job with more cores.

## Input format

The `modisco motifs -i` input is the score file written by `chrombpnet contribs_bw` / the interpret step:
`/raw/seq` (int8 one-hot) and `/shap/seq` (float16 hypothetical contributions), both `(N, 4, L)`, plus
`/projected_shap/seq`, which MoDISco does not read. The datasets are Blosc-compressed; `h5py` needs
`import hdf5plugin` to read them.
