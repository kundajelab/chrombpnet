#!/bin/bash
#SBATCH --job-name=chrombpnet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#
# ChromBPNet 2.x (Keras 3 / JAX) on one GPU: H100 (sm_90), B200 (sm_100), RTX PRO 6000 Blackwell (sm_120).
# Edit the #SBATCH lines for your cluster (partition, gres such as gpu:h100:1, account) and the inputs below.
#
# Once, on a login node (it has no GPU, hence the override; compute nodes then need no internet):
#   git clone https://github.com/kundajelab/chrombpnet && cd chrombpnet
#   CONDA_OVERRIDE_CUDA=13.0 pixi install --locked -e cuda13
# Then, from the checkout:  sbatch workflows/slurm/chrombpnet_gpu.sbatch
#
# NVIDIA driver >= 580 is needed for cuda13 (`nvidia-smi` shows "CUDA Version: 13.x"). With an older driver,
# install and use the cuda12 env instead (PIXI_ENV=cuda12 sbatch ...).
# Do NOT `module load cuda`/`cudnn` and do not set LD_LIBRARY_PATH for this job: JAX ships its own CUDA
# libraries as pip wheels, and system copies found first would shadow them.
set -euo pipefail

REPO="${CHROMBPNET_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
PIXI_ENV="${PIXI_ENV:-cuda13}"

# ---- inputs (edit) ---------------------------------------------------------------------------------------------
FRAGMENTS=/path/to/sample.fragments.tsv.gz      # or use -ibam / -itag below
DATA_TYPE=ATAC                                  # ATAC or DNASE
GENOME=/path/to/hg38.fa
CHROM_SIZES=/path/to/hg38.chrom.sizes
PEAKS=/path/to/peaks.narrowPeak                 # 10 columns
NONPEAKS=/path/to/nonpeaks.bed                  # from `chrombpnet prep nonpeaks`
FOLD=/path/to/fold_0.json
BIAS_MODEL=/path/to/bias.h5                     # chrombpnet 1.x bias models load as-is
OUT=/path/to/output/dir                         # must not exist yet
# -----------------------------------------------------------------------------------------------------------------

if [[ ! -x "${REPO}/.pixi/envs/${PIXI_ENV}/bin/python" ]]; then
    echo "pixi env ${PIXI_ENV} is not installed in ${REPO}. On a login node run:" >&2
    echo "  cd ${REPO} && CONDA_OVERRIDE_CUDA=13.0 pixi install --locked -e ${PIXI_ENV}" >&2
    exit 1
fi
if [[ "${LD_LIBRARY_PATH:-}" == *cuda* || "${LD_LIBRARY_PATH:-}" == *cudnn* ]]; then
    echo "WARNING: LD_LIBRARY_PATH contains CUDA libraries (${LD_LIBRARY_PATH}); they can shadow JAX's own." >&2
fi

nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader
driver_major=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -n1 | cut -d. -f1)
if [[ "${PIXI_ENV}" == cuda13* && "${driver_major}" -lt 580 ]]; then
    echo "NVIDIA driver ${driver_major}.x < 580 cannot run CUDA 13: install and use PIXI_ENV=cuda12" >&2
    exit 1
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"   # TF-MoDISco
export KERAS_BACKEND=jax
# chrombpnet allocates GPU memory on demand. On a GPU shared with other jobs, also cap it:
# export XLA_PYTHON_CLIENT_MEM_FRACTION=0.5
# Caches: keep them off a small $HOME quota if needed.
# export PIXI_CACHE_DIR=/scratch/$USER/pixi-cache JAX_COMPILATION_CACHE_DIR=/scratch/$USER/jax-cache

cd "${REPO}"
# --frozen: use pixi.lock exactly as committed, never rewrite it.
pixi run -e "${PIXI_ENV}" --frozen gpu-check

pixi run -e "${PIXI_ENV}" --frozen chrombpnet pipeline \
    -ifrag "${FRAGMENTS}" \
    -d "${DATA_TYPE}" \
    -g "${GENOME}" \
    -c "${CHROM_SIZES}" \
    -p "${PEAKS}" \
    -n "${NONPEAKS}" \
    -fl "${FOLD}" \
    -b "${BIAS_MODEL}" \
    -o "${OUT}" \
    --device gpu \
    "$@"
# Extra arguments after the script name are passed to chrombpnet, e.g.
#   sbatch workflows/slurm/chrombpnet_gpu.sbatch --optimizer muon
# For a bias model instead: `chrombpnet bias pipeline ... -b 0.5` (bias threshold factor, no -b model).
