# ChromBPNet on a SLURM GPU cluster

`chrombpnet_gpu.sbatch` is a template for one `chrombpnet pipeline` run on one GPU, using the pixi `cuda13`
environment. It targets H100 (sm_90), B200 (sm_100) and RTX PRO 6000 Blackwell (sm_120) nodes. The JAX CUDA 13
wheels ship native code for all three, so there is no PTX JIT warm-up.

## Install once, on a login node

```bash
curl -fsSL https://pixi.sh/install.sh | bash          # if pixi is not available (>= 0.81 recommended)
git clone https://github.com/kundajelab/chrombpnet && cd chrombpnet
CONDA_OVERRIDE_CUDA=13.0 pixi install --locked -e cuda13       # several GB; login nodes have no GPU, hence the override
```

The environment lives in `chrombpnet/.pixi/envs/cuda13`. Compute nodes only read it, so they need no internet.
If `$HOME` has a small quota, set `PIXI_CACHE_DIR` to a scratch directory before installing.

## Driver and CUDA

* `cuda13` needs NVIDIA driver **>= 580**. Check with `nvidia-smi`: the header must say `CUDA Version: 13.x`.
  With an older driver (>= 525), use the `cuda12` environment: `CONDA_OVERRIDE_CUDA=12.0 pixi install --locked -e cuda12`,
  then `PIXI_ENV=cuda12 sbatch ...`. Blackwell GPUs need at least driver 570 even with cuda12.
* **No `module load cuda`/`cudnn` and no `LD_LIBRARY_PATH`.** JAX brings CUDA, cuDNN and NCCL as pip wheels
  inside the environment. System libraries found first through `LD_LIBRARY_PATH` shadow them and cause
  version-mismatch errors. The template warns if `LD_LIBRARY_PATH` mentions CUDA.
* Quick check on a GPU node: `srun --gres=gpu:1 --pty pixi run -e cuda13 gpu-check`, then
  `pixi run -e cuda13-dev test-gpu` for the GPU test suite (install `cuda13-dev` for that).

## Submitting

Edit the `#SBATCH` header (partition, `--gres=gpu:h100:1` or similar, account) and the input paths. Then:

```bash
sbatch workflows/slurm/chrombpnet_gpu.sbatch
sbatch workflows/slurm/chrombpnet_gpu.sbatch --optimizer muon --ema   # extra args go to chrombpnet
```

The script checks that the environment is installed and the driver is new enough. It runs `gpu-check`, then
`chrombpnet pipeline ... --device gpu`, so a job that lands on a node without a working GPU fails immediately
instead of training on the CPU. `pixi run --frozen` uses `pixi.lock` exactly as committed.

## Memory and threads

* JAX allocates GPU memory on demand (`XLA_PYTHON_CLIENT_PREALLOCATE=false`, set by importing chrombpnet when
  unset; `gpu-check` imports chrombpnet first), unlike its default of reserving 75% of the card. It does not return that memory until
  the process ends, and DeepSHAP picks its batch size from the free GPU memory. On a GPU shared with other jobs,
  also set `XLA_PYTHON_CLIENT_MEM_FRACTION` (a fraction of the card's total memory).
* `OMP_NUM_THREADS` and `NUMBA_NUM_THREADS` follow `--cpus-per-task`. TF-MoDISco (numba) is the CPU-heavy step
  at the end of the pipeline, so more CPUs shorten it. `--tomtom-lite` makes the motif report much faster.
* Host memory: the pipeline loads all training sequences, so its memory grows with the number of peaks and
  nonpeaks. After a first run, check its peak memory with `sacct -j <jobid> --format=JobID,MaxRSS,Elapsed` and
  set `--mem` from that.

## Containers instead of pixi

```bash
apptainer pull chrombpnet.sif docker://ghcr.io/kundajelab/chrombpnet:latest-cuda13
apptainer exec --nv --cleanenv --env CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES}" \
    chrombpnet.sif chrombpnet pipeline ...
```

The images come from the repository's docker workflow (`-cuda13` and `-cuda12` tags) once the maintainers enable
it. Until then, build one with `docker build` (see the `Dockerfile` header) and convert it with
`apptainer build chrombpnet.sif docker-daemon://chrombpnet:cuda13`.

`--nv` injects the host driver. The image needs no CUDA module and sets `KERAS_BACKEND=jax` itself. It also
ignores a host `PYTHONPATH`/`PYTHONHOME` (e.g. from `module load python`) and the packages in `~/.local`, which
Apptainer would otherwise see through the bound `$HOME`. `--cleanenv` keeps out the rest of the job's
environment, including the `CUDA_VISIBLE_DEVICES` that SLURM sets for the allocated GPUs, hence the `--env`.
Drop that `--env` where the scheduler does not set the variable (without `--gres=gpu`): an empty
`CUDA_VISIBLE_DEVICES` hides every GPU. Forward other variables the same way, e.g.
`--env XLA_PYTHON_CLIENT_MEM_FRACTION=0.5` on a shared GPU.
