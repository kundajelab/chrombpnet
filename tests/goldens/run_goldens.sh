#!/usr/bin/env bash
# Produce legacy parity goldens, CPU only, inside an apptainer container with chrombpnet 1.x (TF-Keras 2.x).
# Read-only on the container and the results tree; writes only to $OUT.
#
# Required environment:
#   SANDBOX  apptainer image or sandbox directory with chrombpnet 1.x, TF 2.x and kundajelab-shap installed
#   RESULTS  chrombpnet 1.x results tree (default: $CHROMBPNET_LEGACY_RESULTS), laid out as
#              bias_models/bias_model_<b>/<PREFIX>/{models,evaluation,auxiliary}/   one directory per bias model <b>
#              full_models/<PREFIX>/{models,auxiliary}/
#   PREFIX   run name (default: $CHROMBPNET_LEGACY_PREFIX): the run directory above, and the file prefix of the
#            bias-model files (<PREFIX>_bias.h5, <PREFIX>_bias_predictions.h5, ...)
#   BIAS     the bias model <b> whose regions give the prediction inputs and whose DeepSHAP goldens are written
#            (default: $CHROMBPNET_LEGACY_BIAS); every bias_models/bias_model_*/ gets prediction goldens
#   GENOME   reference FASTA (with its .fai) the models were trained on
# Optional: OUT (default ./goldens), BIND (comma-separated paths to bind into the container; default: the
# directories of RESULTS, GENOME, OUT and this script). Runs on $SLURM_CPUS_PER_TASK threads (default 2).
# Needs `apptainer` on PATH; for example:
#   SANDBOX=<chrombpnet 1.x .sif or sandbox dir> RESULTS=<results tree> PREFIX=<run name> BIAS=<b> \
#       GENOME=<genome.fa> OUT=<scratch dir>/goldens bash <checkout>/tests/goldens/run_goldens.sh
set -euo pipefail

RESULTS=${RESULTS:-${CHROMBPNET_LEGACY_RESULTS:-}}
PREFIX=${PREFIX:-${CHROMBPNET_LEGACY_PREFIX:-}}
BIAS=${BIAS:-${CHROMBPNET_LEGACY_BIAS:-}}
missing=()
for v in SANDBOX RESULTS PREFIX BIAS GENOME; do
    [[ -n ${!v:-} ]] || missing+=("$v")
done
if (( ${#missing[@]} )); then
    echo "run_goldens.sh: set ${missing[*]} (see the header of $0)" >&2
    exit 2
fi
die() { echo "run_goldens.sh: $*" >&2; exit 2; }

R=$RESULTS
OUT=${OUT:-$PWD/goldens}
SCRIPT=$(cd "$(dirname "$0")" && pwd)/make_goldens.py
BM=$R/bias_models
FM=$R/full_models/$PREFIX
BR=$BM/bias_model_$BIAS/$PREFIX
BIAS_H5=$BR/models/${PREFIX}_bias.h5
NT=${SLURM_CPUS_PER_TASK:-2}
[[ -d $R ]] || die "RESULTS=$R is not a directory"
[[ -f $GENOME ]] || die "GENOME=$GENOME is not a file"
[[ -d $FM/models ]] || die "no full_models/$PREFIX/models under RESULTS=$R (check PREFIX)"
[[ -f $BIAS_H5 ]] || die "no $BIAS_H5 (check BIAS and PREFIX)"
mkdir -p "$OUT" "$OUT/.mpl"
BIND=${BIND:-$R,$(dirname "$GENOME"),$OUT,$(dirname "$SCRIPT")}

export APPTAINERENV_CUDA_VISIBLE_DEVICES=""
export APPTAINERENV_TF_CPP_MIN_LOG_LEVEL=2
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
export APPTAINERENV_PYTHONNOUSERSITE=1
export APPTAINERENV_PYTHONUNBUFFERED=1
export APPTAINERENV_MPLCONFIGDIR=$OUT/.mpl
export APPTAINERENV_OMP_NUM_THREADS=$NT
export APPTAINERENV_TF_NUM_INTRAOP_THREADS=$NT
export APPTAINERENV_TF_NUM_INTEROP_THREADS=1
unset PYTHONPATH PYTHONHOME VIRTUAL_ENV

legacy() { apptainer exec --cleanenv --bind "$BIND" "$SANDBOX" python3 "$SCRIPT" "$@"; }

legacy_versions() {
    apptainer exec --cleanenv --bind "$BIND" "$SANDBOX" python3 -c '
import json, sys, numpy, tensorflow as tf, shap, chrombpnet
print(json.dumps({"python": sys.version.split()[0], "tensorflow": tf.__version__, "numpy": numpy.__version__,
                  "shap": getattr(shap, "__version__", "?"), "chrombpnet": chrombpnet.__file__}))'
}

START=$(date -u +%FT%TZ)
legacy_versions > "$OUT/versions.json"
cat "$OUT/versions.json"

bias_models=()
for d in "$BM"/bias_model_*/; do
    b=${d%/}
    b=${b##*/bias_model_}
    f=$BM/bias_model_$b/$PREFIX/models/${PREFIX}_bias.h5
    if [[ -f $f ]]; then bias_models+=("bias_$b=$f"); fi
done

echo "== predict"
legacy predict --out "$OUT" --genome "$GENOME" \
    --coords-h5 "$BR/evaluation/${PREFIX}_bias_predictions.h5" \
    --models "${bias_models[@]}" "bias_model_scaled=$FM/models/bias_model_scaled.h5" \
             "chrombpnet=$FM/models/chrombpnet.h5" "chrombpnet_nobias=$FM/models/chrombpnet_nobias.h5"

echo "== shap bias_$BIAS (counts + profile)"
legacy shap --out "$OUT" --genome "$GENOME" --name "bias_$BIAS" --model "$BIAS_H5" \
    --regions "$BR/auxiliary/interpret_subsample/${PREFIX}_bias.interpreted_regions.bed" \
    --n 64 --heads counts profile
legacy fvals --out "$OUT" --name "bias_$BIAS" --model "$BIAS_H5"

echo "== shap chrombpnet_nobias (profile)"
legacy shap --out "$OUT" --genome "$GENOME" --name chrombpnet_nobias --model "$FM/models/chrombpnet_nobias.h5" \
    --regions "$FM/auxiliary/30K_subsample_peaks.bed" --n 16 --heads profile counts
legacy fvals --out "$OUT" --name chrombpnet_nobias --model "$FM/models/chrombpnet_nobias.h5"

echo "== traces"
legacy trace --out "$OUT" --name bias_128x4 --bigwig "$FM/auxiliary/data_unstranded.bw" \
    --filters 128 --n-dil-layers 4 --counts-loss-weight 10
legacy trace --out "$OUT" --name chrombpnet_32x4 --bigwig "$FM/auxiliary/data_unstranded.bw" \
    --filters 32 --n-dil-layers 4 --counts-loss-weight 10 --bias-model "$FM/models/bias_model_scaled.h5"

python3 - "$OUT" "$START" <<'PY'
import hashlib, json, os, sys
out, start = sys.argv[1], sys.argv[2]
files = {}
for root, _, names in os.walk(out):
    for n in names:
        p = os.path.join(root, n)
        if "/.mpl" in p or n == "manifest.json":
            continue
        h = hashlib.md5(open(p, "rb").read()).hexdigest()
        files[os.path.relpath(p, out)] = {"bytes": os.path.getsize(p), "md5": h}
json.dump({"started": start, "versions": json.load(open(os.path.join(out, "versions.json"))), "files": files},
          open(os.path.join(out, "manifest.json"), "w"), indent=1)
print("wrote manifest with", len(files), "files")
PY
du -sh "$OUT"
