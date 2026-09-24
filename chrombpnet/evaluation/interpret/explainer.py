"""Native DeepSHAP explainer for ChromBPNet / BPNet Keras 3 models on the JAX backend.

Usage::

    explainer = DeepLiftShap(model, heads=("counts", "profile"))
    hyp = explainer.explain(onehot_seqs)          # {"counts": (N, L, 4), "profile": (N, L, 4)} float32

Multipliers are the kundajelab-shap TFDeepExplainer ones (see deeplift_jax): ReLU activations of Conv1D/Dense
layers are swapped for the rescale rule while the model is traced, one jitted forward per chunk of S sequences
x K references is shared by the requested heads, and the hypothetical projection runs on the device.
"""
import contextlib
import time
import warnings

import jax
import jax.numpy as jnp
import keras
import numpy as np

from chrombpnet.evaluation.interpret import deeplift_jax
from chrombpnet.evaluation.interpret.dinuc_shuffle import NUM_SHUFFLES, make_references
from chrombpnet.training.utils.layers import LogSumExp

HEADS = ("counts", "profile")
PRECISIONS = ("auto", "highest", "high", "default")
ADDITIVITY_THRESHOLD = 1e-3

_LINEAR_LAYERS = (
    keras.layers.InputLayer,
    keras.layers.Cropping1D,
    keras.layers.Add,
    keras.layers.Flatten,
    keras.layers.GlobalAveragePooling1D,
    keras.layers.Concatenate,
)
_ACTIVATION_LAYERS = (keras.layers.Conv1D, keras.layers.Dense)


def iter_layers(model):
    """All layers of a model, recursing into nested models (e.g. the two submodels of a full chrombpnet.h5)."""
    for layer in model.layers:
        yield layer
        if isinstance(layer, keras.Model):
            yield from iter_layers(layer)


def activation_name(layer):
    act = getattr(layer, "activation", None)
    if act is None:
        return "linear"
    try:
        name = keras.activations.serialize(act)
    except Exception:
        name = getattr(act, "__name__", None)
    return name if isinstance(name, str) else repr(name)


def check_model(model, heads):
    """Raise NotImplementedError unless every layer has a DeepLIFT rule here (ReLU rescale, or linear)."""
    if len(model.inputs) != 1:
        raise NotImplementedError("DeepSHAP supports single-input models only ({} has {} inputs)".format(
            model.name, len(model.inputs)))
    n_out = len(model.outputs)
    for head in heads:
        if head not in HEADS:
            raise ValueError("unknown head {!r}; expected one of {}".format(head, HEADS))
        if head == "counts" and n_out < 2:
            raise ValueError("counts head requested but the model has a single output")
    has_logsumexp = False
    for layer in iter_layers(model):
        if isinstance(layer, keras.Model):
            continue
        if isinstance(layer, LogSumExp):
            has_logsumexp = True
        elif isinstance(layer, _ACTIVATION_LAYERS):
            act = activation_name(layer)
            if act not in ("relu", "linear"):
                raise NotImplementedError(
                    "layer {!r} ({}) uses activation {!r}; DeepSHAP here only has rules for relu and linear".format(
                        layer.name, type(layer).__name__, act))
        elif not isinstance(layer, _LINEAR_LAYERS):
            raise NotImplementedError(
                "layer {!r} ({}) has no DeepLIFT rule; supported layers: InputLayer, Conv1D and Dense (relu or "
                "linear), Cropping1D, Add, Concatenate, Flatten, GlobalAveragePooling1D, nested models".format(
                    layer.name, type(layer).__name__))
    if has_logsumexp and "counts" in heads:
        raise NotImplementedError(
            "counts-head DeepSHAP of a full chrombpnet model (logsumexp of the bias and no-bias counts) is not "
            "supported: interpret the bias-corrected *_nobias.h5 model instead (as the chrombpnet pipeline does). "
            "The profile head of the full model is supported.")


@contextlib.contextmanager
def deeplift_rules(model, relu=deeplift_jax.dl_relu):
    """Temporarily give every relu Conv1D/Dense the DeepLIFT rule, and run every layer in float32.

    Enter it while the model is traced (inside the jitted function) so any retrace also sees the rules.
    """
    swapped, policies = [], []
    try:
        for layer in [model] + list(iter_layers(model)):
            if isinstance(layer, _ACTIVATION_LAYERS) and activation_name(layer) == "relu":
                swapped.append((layer, layer.activation))
                layer.activation = relu
            if layer.compute_dtype in ("bfloat16", "float16"):
                policies.append((layer, layer.dtype_policy))
                layer.dtype_policy = "float32"
        yield
    finally:
        for layer, act in reversed(swapped):
            layer.activation = act
        for layer, policy in reversed(policies):
            layer.dtype_policy = policy


def _activation_floats_per_row(model):
    """Conv1D output values per input row, summed over layers (what the backward pass keeps per row)."""
    total = 0
    for layer in iter_layers(model):
        if isinstance(layer, keras.layers.Conv1D):
            shape = layer.output.shape[1:]
            if any(d is None for d in shape):
                return None
            total += int(np.prod(shape))
    return total


def auto_batch_seqs(model, n_heads=1, num_shuffles=NUM_SHUFFLES, device=None):
    """Sequences per step, sized from the model's activations and the free device memory.

    Each sequence contributes 2 * num_shuffles rows (itself repeated and its references). On an RTX PRO 6000
    (TF32), throughput stops improving once a step holds a few sequences, for chrombpnet_nobias (512x8) as for a
    128x4 bias model, so steps are kept to ~12 GB (`budget` below) with a conservative per-row estimate.
    """
    per_row = None
    try:
        per_row = _activation_floats_per_row(model)
    except Exception:  # symbolic shapes unavailable (e.g. a subclassed model)
        pass
    if not per_row:
        filters = [layer.filters for layer in iter_layers(model) if isinstance(layer, keras.layers.Conv1D)]
        return 4 if filters and max(filters) >= 256 else 16
    bytes_per_seq = 2 * num_shuffles * per_row * 4 * (3.0 + 1.5 * (n_heads - 1))
    budget = 12e9
    device = device or jax.devices()[0]
    try:
        stats = device.memory_stats() or {}
    except Exception:
        stats = {}
    if stats.get("bytes_limit"):
        free = stats["bytes_limit"] - stats.get("bytes_in_use", 0)
        budget = min(budget, 0.6 * free)
    return int(max(1, min(32, budget // bytes_per_seq)))


def resolve_precision(precision, backend=None):
    """'auto' -> full float32 on CPU (free there) and the backend default (TF32 on Ampere+) on GPU.

    Full float32 convolutions are extremely slow on some GPUs (far slower than TF32 on an RTX PRO 6000
    Blackwell), and chrombpnet 1.x also ran DeepSHAP in TF32 on GPUs.
    """
    if precision != "auto":
        return precision
    return "highest" if (backend or jax.default_backend()) == "cpu" else "default"


def precision_context(precision):
    if precision not in PRECISIONS:
        raise ValueError("precision must be one of {} for DeepSHAP, got {!r}".format(PRECISIONS, precision))
    if precision == "default":
        return contextlib.nullcontext()
    return jax.default_matmul_precision(precision)


def _is_oom(err):
    msg = str(err)
    return "RESOURCE_EXHAUSTED" in msg or "Out of memory" in msg or "out of memory" in msg


class ExplainResult:
    """hyp[head]: (N, L, 4) float32; mult[head]: (N, K, L, 4) float32 (if requested); refs: (N, K, L, 4) int8
    (if requested); additivity[head]: summary of the per-pair summation-to-delta check."""

    def __init__(self):
        self.hyp, self.mult, self.additivity = {}, {}, {}
        self.refs = None


class DeepLiftShap:
    """DeepSHAP (DeepLIFT rescale rule averaged over dinucleotide-shuffled references) for BPNet-style models.

    Args:
        model: Keras 3 model with outputs [profile logits, log counts].
        heads: subset of ("counts", "profile").
        profile_weighting: see deeplift_jax.profile_weights; "chrombpnet" reproduces chrombpnet 1.x.
        precision: matmul/conv precision while tracing: "auto" (default: "highest" on CPU, "default" on GPU),
            "highest" (full float32), "high", or "default" (TF32 on Ampere+ GPUs; faster, noisier near ReLU
            branch switches).
        batch_seqs: sequences per step (each with num_shuffles references); None picks auto_batch_seqs(model).
        num_shuffles: references per sequence when they are generated here.
        additivity_threshold: warn when a pair's |sum m*(x-r) - delta| / max(|delta|, 1) exceeds it.
        variant: "paired" (reference-half cotangent zeroed) or "full" (TFDeepExplainer-like, for tests).
    """

    def __init__(self, model, heads=HEADS, profile_weighting="chrombpnet", precision="auto", batch_seqs=None,
                 num_shuffles=NUM_SHUFFLES, additivity_threshold=ADDITIVITY_THRESHOLD, variant="paired"):
        heads = tuple(h for h in HEADS if h in tuple(heads))
        if not heads:
            raise ValueError("no heads requested")
        check_model(model, heads)
        if profile_weighting not in deeplift_jax.PROFILE_WEIGHTINGS:
            raise ValueError("profile_weighting must be one of {}, got {!r}".format(
                deeplift_jax.PROFILE_WEIGHTINGS, profile_weighting))
        precision = resolve_precision(precision)
        precision_context(precision)
        if variant not in ("paired", "full"):
            raise ValueError("variant must be 'paired' or 'full'")
        self.model = model
        self.heads = heads
        self.profile_weighting = profile_weighting
        self.precision = precision
        self.batch_seqs = int(batch_seqs) if batch_seqs else None
        self.num_shuffles = int(num_shuffles)
        self.additivity_threshold = additivity_threshold
        self.variant = variant
        self.input_length = model.inputs[0].shape[1]
        self._steps = {}

    # ------------------------------------------------------------------ jitted step
    def _step(self, n_seqs, n_refs, seq_len, return_mult):
        key = (n_seqs, n_refs, seq_len, return_mult)
        if key not in self._steps:
            self._steps[key] = jax.jit(self._make_step(n_seqs, n_refs, seq_len, return_mult))
        return self._steps[key]

    def _make_step(self, n_seqs, n_refs, seq_len, return_mult):
        model, heads, mode = self.model, self.heads, self.profile_weighting
        full = self.variant == "full"
        relu = deeplift_jax.dl_relu_full if full else deeplift_jax.dl_relu
        n_pairs = n_seqs * n_refs

        def forward(tv, ntv, joint):
            with deeplift_rules(model, relu):
                outputs, _ = model.stateless_call(tv, ntv, joint, training=False)
            if not isinstance(outputs, (list, tuple)):
                outputs = [outputs]
            logits = outputs[0]
            logcounts = outputs[1] if len(outputs) > 1 else jnp.zeros((joint.shape[0], 1), logits.dtype)
            return deeplift_jax.targets(logits, logcounts, mode=mode, full=full)

        def step(tv, ntv, x, refs):
            xf = x.astype(jnp.float32)
            rf = refs.astype(jnp.float32)
            x_rows = jnp.repeat(xf, n_refs, axis=0)
            r_rows = jnp.reshape(rf, (n_pairs, seq_len, xf.shape[-1]))
            joint = jnp.concatenate([x_rows, r_rows], axis=0)
            t, vjp_fn, deltas = jax.vjp(lambda j: forward(tv, ntv, j), joint, has_aux=True)
            out = {}
            for head in heads:
                cot = (jnp.ones_like(t[0]), jnp.zeros_like(t[1])) if head == "counts" else \
                    (jnp.zeros_like(t[0]), jnp.ones_like(t[1]))
                m = vjp_fn(cot)[0][:n_pairs].astype(jnp.float32)
                m4 = jnp.reshape(m, (n_seqs, n_refs, seq_len, -1))
                res = {
                    "hyp": deeplift_jax.hypothetical(m4, rf),
                    "sum": jnp.reshape(deeplift_jax.summation_delta(m, x_rows, r_rows), (n_seqs, n_refs)),
                    "delta": jnp.reshape(deltas[head].astype(jnp.float32), (n_seqs, n_refs)),
                }
                if return_mult:
                    res["mult"] = m4
                out[head] = res
            return out

        return step

    # ------------------------------------------------------------------ driver
    def _reference_fn(self, references, seed):
        if references is None:
            return lambda x, start, stop: make_references(x, num_shuffles=self.num_shuffles, seed=seed)
        if callable(references):
            return references
        refs = np.asarray(references)
        return lambda x, start, stop: refs[start:stop]

    def iter_explain(self, seqs, references=None, seed=1234, return_multipliers=False, additivity=None):
        """Yield (start, stop, refs, {head: {"hyp", ["mult"]}}) chunk by chunk, in input order.

        references: None (num_shuffles content-seeded dinucleotide shuffles per sequence, see
        dinuc_shuffle.make_references), an (N, K, L, 4) array, or a callable (x_chunk, start, stop) -> refs.
        additivity: optional dict, filled with the per-head additivity summary.
        """
        seqs = np.asarray(seqs)
        if seqs.ndim != 3 or seqs.shape[2] != 4:
            raise ValueError("expected one-hot sequences of shape (N, L, 4), got {}".format(seqs.shape))
        if self.input_length is not None and seqs.shape[1] != self.input_length:
            raise ValueError("sequences have length {} but the model expects {}".format(
                seqs.shape[1], self.input_length))
        n = seqs.shape[0]
        stats = additivity if additivity is not None else {}
        for head in self.heads:
            stats[head] = {"max_rel_err": 0.0, "max_abs_err": 0.0, "n_pairs": 0, "n_over_threshold": 0,
                           "threshold": self.additivity_threshold}
        if n == 0:
            return
        get_refs = self._reference_fn(references, seed)
        tv = [v.value for v in self.model.trainable_variables]
        ntv = [v.value for v in self.model.non_trainable_variables]
        chunk = min(self.batch_seqs or auto_batch_seqs(self.model, len(self.heads), self.num_shuffles), n)

        def load(start):
            stop = min(start + chunk, n)
            x = np.ascontiguousarray(seqs[start:stop], dtype=np.int8)
            refs = np.ascontiguousarray(get_refs(x, start, stop), dtype=np.int8)
            if refs.ndim != 4 or refs.shape[0] != stop - start or refs.shape[2:] != x.shape[1:]:
                raise ValueError("references for sequences {}:{} have shape {}, expected ({}, K, {}, {})".format(
                    start, stop, refs.shape, stop - start, x.shape[1], x.shape[2]))
            return start, stop, x, refs

        start = 0
        nxt = load(0)
        while start < n:
            if nxt is None or nxt[0] != start or nxt[1] != min(start + chunk, n):
                nxt = load(start)
            _, stop, x, refs = nxt
            n_refs = refs.shape[1]
            pad = chunk - (stop - start)
            if pad:
                x_in = np.concatenate([x, np.repeat(x[-1:], pad, axis=0)])
                r_in = np.concatenate([refs, np.repeat(refs[-1:], pad, axis=0)])
            else:
                x_in, r_in = x, refs
            try:
                # the precision is read while tracing and is part of the jit cache key: same context every call
                with precision_context(self.precision):
                    out = self._step(chunk, n_refs, x.shape[1], return_multipliers)(tv, ntv, x_in, r_in)
                # prepare the next chunk's references on the host while the device works
                nxt = load(stop) if stop < n else None
                out = jax.device_get(out)
            except Exception as err:  # XlaRuntimeError surfaces at dispatch or at device_get
                if not _is_oom(err) or chunk == 1:
                    raise
                chunk = max(1, chunk // 2)
                warnings.warn("DeepSHAP ran out of device memory; retrying with {} sequences per step".format(chunk))
                nxt = None
                continue
            valid = stop - start
            results = {}
            for head in self.heads:
                res = out[head]
                self._update_additivity(stats[head], res["sum"][:valid], res["delta"][:valid])
                results[head] = {"hyp": np.asarray(res["hyp"][:valid])}
                if return_multipliers:
                    results[head]["mult"] = np.asarray(res["mult"][:valid])
            yield start, stop, refs, results
            start = stop
        self._report_additivity(stats)

    @staticmethod
    def _update_additivity(stat, total, delta):
        total = np.asarray(total, dtype=np.float64)
        delta = np.asarray(delta, dtype=np.float64)
        abs_err = np.abs(total - delta)
        rel_err = abs_err / np.maximum(np.abs(delta), 1.0)
        stat["n_pairs"] += int(rel_err.size)
        if rel_err.size:
            stat["max_rel_err"] = max(stat["max_rel_err"], float(rel_err.max()))
            stat["max_abs_err"] = max(stat["max_abs_err"], float(abs_err.max()))
            stat["n_over_threshold"] += int((rel_err > stat["threshold"]).sum())

    def _report_additivity(self, stats):
        for head, stat in stats.items():
            if stat["n_over_threshold"]:
                warnings.warn(
                    "DeepSHAP {} head: {} of {} (input, reference) pairs violate summation-to-delta by more than "
                    "{:g} (max relative error {:.3g}); expected only with precision='default' on GPUs (TF32) or "
                    "for layers without a DeepLIFT rule".format(head, stat["n_over_threshold"], stat["n_pairs"],
                                                                 stat["threshold"], stat["max_rel_err"]))

    def run(self, seqs, references=None, seed=1234, return_multipliers=False, return_references=False):
        """Explain all sequences at once; returns an ExplainResult."""
        seqs = np.asarray(seqs)
        result = ExplainResult()
        parts = {head: [] for head in self.heads}
        mults = {head: [] for head in self.heads}
        refs_parts = []
        for _, _, refs, res in self.iter_explain(seqs, references, seed, return_multipliers, result.additivity):
            for head in self.heads:
                parts[head].append(res[head]["hyp"])
                if return_multipliers:
                    mults[head].append(res[head]["mult"])
            if return_references:
                refs_parts.append(refs)
        for head in self.heads:
            result.hyp[head] = np.concatenate(parts[head]) if parts[head] else \
                np.zeros((0,) + seqs.shape[1:], np.float32)
            if return_multipliers:
                result.mult[head] = np.concatenate(mults[head]) if mults[head] else \
                    np.zeros((0, self.num_shuffles) + seqs.shape[1:], np.float32)
        if return_references:
            result.refs = np.concatenate(refs_parts) if refs_parts else \
                np.zeros((0, self.num_shuffles) + seqs.shape[1:], np.int8)
        return result

    def explain(self, seqs, references=None, seed=1234):
        """Hypothetical contributions {head: (N, L, 4) float32}; multiply by the one-hot input to project."""
        return self.run(seqs, references, seed).hyp

    def multipliers(self, seqs, references=None, seed=1234):
        """Raw input-half multipliers {head: (N, K, L, 4) float32} (the TFDeepExplainer `mult`)."""
        return self.run(seqs, references, seed, return_multipliers=True).mult


def timed_iter(iterator, total, every=60.0, label="DeepSHAP"):
    """Pass chunks through, printing progress about every `every` seconds."""
    t0 = last = time.time()
    for item in iterator:
        yield item
        now = time.time()
        done = item[1]
        if now - last >= every or done == total:
            rate = done / max(now - t0, 1e-9)
            print("{}: {}/{} sequences ({:.1f} seq/s)".format(label, done, total, rate), flush=True)
            last = now
