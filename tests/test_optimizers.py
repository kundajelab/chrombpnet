"""make_optimizer (Adam / ChromBPNetMuon, EMA, cosine schedule) and the runtime precision/device helpers."""
import json
import math
import os
import subprocess
import sys
import types

import numpy as np
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import jax
import keras

import chrombpnet.training.models.bpnet_model as bpnet_model
import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
from chrombpnet.training import optimizers, runtime
from chrombpnet.training.optimizers import ChromBPNetMuon, make_optimizer
from chrombpnet.training.utils import model_io

INPUTLEN, OUTPUTLEN = 514, 200


def params(n_dil_layers=3, **kw):
    p = {"filters": "8", "n_dil_layers": str(n_dil_layers), "counts_loss_weight": "10", "inputlen": str(INPUTLEN),
         "outputlen": str(OUTPUTLEN)}
    p.update(kw)
    return p


def muon_routed(optimizer, variables):
    """Paths of the variables that get the Muon update (the ones without an Adam velocity slot)."""
    optimizer.build(variables)
    return [v.path for v in variables if optimizer.adam_velocities[optimizer._get_variable_index(v)] is None]


def test_default_is_legacy_adam():
    opt = make_optimizer(types.SimpleNamespace())
    assert type(opt) is keras.optimizers.Adam
    assert float(opt.learning_rate) == pytest.approx(1e-3)
    assert (opt.beta_1, opt.beta_2, opt.epsilon, opt.amsgrad) == (0.9, 0.999, 1e-7, False)
    assert not opt.use_ema and opt.weight_decay is None
    opt = make_optimizer(types.SimpleNamespace(learning_rate=5e-4, optimizer=None, ema=None, lr_schedule=None))
    assert type(opt) is keras.optimizers.Adam and float(opt.learning_rate) == pytest.approx(5e-4)


def test_ema_options():
    opt = make_optimizer(types.SimpleNamespace(ema=True))
    assert opt.use_ema and opt.ema_momentum == pytest.approx(0.999)
    opt = make_optimizer(types.SimpleNamespace(optimizer="muon", ema=True, ema_momentum=0.99))
    assert isinstance(opt, ChromBPNetMuon) and opt.use_ema and opt.ema_momentum == pytest.approx(0.99)


def test_muon_settings():
    opt = make_optimizer(types.SimpleNamespace(optimizer="muon", learning_rate=1e-3))
    assert isinstance(opt, ChromBPNetMuon)
    assert float(opt.learning_rate) == pytest.approx(optimizers.DEFAULT_MUON_LR)
    assert float(opt.learning_rate) * opt.adam_lr_ratio == pytest.approx(1e-3)
    assert opt.weight_decay is None and opt.adam_weight_decay is None
    assert (opt.momentum, opt.nesterov, opt.ns_steps, opt.rms_rate) == (0.95, True, 5, 0.2)
    opt = make_optimizer(types.SimpleNamespace(optimizer="muon", learning_rate=1e-3, muon_lr=4e-3))
    assert float(opt.learning_rate) == pytest.approx(4e-3) and opt.adam_lr_ratio == pytest.approx(0.25)
    restored = keras.optimizers.deserialize(keras.optimizers.serialize(opt))
    assert isinstance(restored, ChromBPNetMuon) and restored.muon_variables == opt.muon_variables


@pytest.mark.parametrize("bad", [dict(optimizer="sgd"), dict(lr_schedule="step"), dict(lr_schedule="cosine"),
                                 dict(lr_schedule="cosine", total_steps=10, warmup_steps=10)])
def test_invalid_options(bad):
    with pytest.raises(ValueError):
        make_optimizer(types.SimpleNamespace(**bad))


def test_cosine_schedule():
    peak, total = 1e-3, 100
    opt = make_optimizer(types.SimpleNamespace(lr_schedule="cosine", total_steps=total))
    schedule = opt._learning_rate
    lr = lambda step: float(schedule(step))
    assert lr(0) == pytest.approx(0.01 * peak)
    assert lr(5) == pytest.approx(0.01 * peak + 0.5 * 0.99 * peak)
    assert lr(10) == pytest.approx(peak)
    assert lr(55) == pytest.approx(peak * (0.01 + 0.99 * 0.5 * (1 + math.cos(math.pi * 0.5))))
    assert lr(total) == pytest.approx(0.01 * peak) and lr(10 * total) == pytest.approx(0.01 * peak)
    assert all(lr(s) >= lr(s + 1) for s in range(10, total))
    opt = make_optimizer(types.SimpleNamespace(lr_schedule="cosine", total_steps=total, warmup_steps=0))
    assert float(opt._learning_rate(0)) == pytest.approx(peak)
    # muon: the schedule drives the Muon lr, the Adam-routed variables follow it through adam_lr_ratio
    opt = make_optimizer(types.SimpleNamespace(optimizer="muon", lr_schedule="cosine", total_steps=total))
    assert float(opt._learning_rate(10)) * opt.adam_lr_ratio == pytest.approx(peak)


def test_muon_routes_only_dilated_conv_kernels(tmp_path):
    bias = bpnet_model.getModelGivenModelOptionsAndWeightInits(
        types.SimpleNamespace(seed=1, learning_rate=1e-3, optimizer="muon"), params(n_dil_layers=3))
    assert isinstance(bias.optimizer, ChromBPNetMuon)
    assert muon_routed(bias.optimizer, bias.trainable_variables) == [
        "bpnet_1conv/kernel", "bpnet_2conv/kernel", "bpnet_3conv/kernel"]

    bias.save(str(tmp_path / "bias.h5"))
    model = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
        types.SimpleNamespace(seed=1, learning_rate=1e-3, optimizer="muon"),
        params(n_dil_layers=2, bias_model_path=str(tmp_path / "bias.h5")))
    # the frozen bias model never reaches the optimizer; first conv, heads and biases stay on Adam
    assert muon_routed(model.optimizer, model.trainable_variables) == [
        "wo_bias_bpnet_1conv/kernel", "wo_bias_bpnet_2conv/kernel"]

    # the optimizer config is written into the .h5; loading never needs it (compile=False)
    model.save(str(tmp_path / "muon.h5"))
    assert model_io.load_model(str(tmp_path / "muon.h5")).output_names == model.output_names


def newton_schulz5(g, steps=5, a=3.4445, b=-4.7750, c=2.0315):
    x = g.astype(np.float64)
    transposed = x.shape[0] > x.shape[1]
    if transposed:
        x = x.T
    x = x / (np.linalg.norm(x) + 1e-7)
    for _ in range(steps):
        s = x @ x.T
        x = a * x + (b * s + c * s @ s) @ x
    return x.T if transposed else x


def test_muon_update_math():
    rng = np.random.RandomState(0)
    dilated = keras.layers.Conv1D(6, 3, name="bpnet_1conv")
    dilated.build((None, 32, 5))
    first = keras.layers.Conv1D(6, 3, name="bpnet_1st_conv")
    first.build((None, 32, 5))
    variables = [dilated.kernel, dilated.bias, first.kernel, first.bias]
    start = [np.array(v) for v in variables]

    args = types.SimpleNamespace(optimizer="muon", learning_rate=1e-3, muon_lr=2e-3)
    muon = make_optimizer(args)
    assert muon_routed(muon, variables) == ["bpnet_1conv/kernel"]
    adam_vars = [keras.Variable(w) for w in start[1:]]
    adam = keras.optimizers.Adam(learning_rate=1e-3)
    grads = [[rng.normal(size=w.shape).astype(np.float32) for w in start] for _ in range(3)]
    for step, g in enumerate(grads):
        muon.apply_gradients(zip(g, variables))
        adam.apply_gradients(zip(g[1:], adam_vars))
        if step == 0:
            # first step: momentum = g, nesterov g + 0.95 g, NS on the (k*in, out) flattening, rms_rate scaling
            flat = (1.95 * g[0]).reshape(-1, g[0].shape[-1])
            expected = start[0] - 2e-3 * 0.2 * math.sqrt(max(flat.shape)) * newton_schulz5(flat).reshape(g[0].shape)
            np.testing.assert_allclose(np.array(dilated.kernel), expected, rtol=2e-4, atol=1e-6)
    # every other variable follows plain Adam at args.learning_rate exactly
    for v, ref in zip(variables[1:], adam_vars):
        np.testing.assert_allclose(np.array(v), np.array(ref), rtol=1e-6, atol=1e-7)


def test_muon_reduces_loss_on_toy_problem():
    args = types.SimpleNamespace(seed=2, learning_rate=1e-3, optimizer="muon", muon_lr=2e-3)
    model = bpnet_model.getModelGivenModelOptionsAndWeightInits(args, params(n_dil_layers=2))
    rng = np.random.RandomState(0)
    x = np.eye(4, dtype=np.float32)[rng.randint(0, 4, (16, INPUTLEN))]
    # profile and counts that depend on the sequence (GC content), so there is something to learn
    y = rng.poisson(1 + 4 * x[:, 157:357, 1:3].sum(-1)).astype(np.float32)
    targets = (y, np.log1p(y.sum(-1, keepdims=True)))
    losses = []
    for _ in range(40):
        model.reset_metrics()
        losses.append(model.train_on_batch(x, targets, return_dict=True)["loss"])
    assert np.isfinite(losses).all()
    assert np.mean(losses[-5:]) < 0.8 * np.mean(losses[:5])


def test_configure_precision():
    try:
        assert runtime.configure_precision("highest") == "highest"
        assert jax.config.jax_default_matmul_precision == "highest"
        # 'default' resets a 'highest' left by an earlier run in the same process
        assert runtime.configure_precision("default") == "default"
        assert jax.config.jax_default_matmul_precision is None
        runtime.configure_precision("highest")
        assert runtime.configure_precision(None) == "default"
        assert jax.config.jax_default_matmul_precision is None
    finally:
        jax.config.update("jax_default_matmul_precision", None)
    previous = keras.config.dtype_policy()
    try:
        runtime.configure_precision("bf16")
        assert keras.config.dtype_policy().name == "mixed_bfloat16"
        assert runtime.bf16_active() and runtime.head_dtype() == "float32"
    finally:
        keras.config.set_dtype_policy(previous)
    assert runtime.configure_precision(None) == "default" and runtime.head_dtype() is None
    assert jax.config.jax_default_matmul_precision is None
    with pytest.raises(ValueError):
        runtime.configure_precision("fp16")


def test_assert_gpu_if_requested():
    jax.devices()  # the backends exist in this process: 'auto' reports them
    backend = runtime.assert_gpu_if_requested("auto")
    assert backend == jax.default_backend()
    assert runtime.assert_gpu_if_requested("auto", initialize=True) == backend
    with pytest.raises(ValueError):
        runtime.assert_gpu_if_requested("tpu")
    if jax.default_backend() == "cpu":
        with pytest.raises(RuntimeError, match="--device gpu"):
            runtime.assert_gpu_if_requested("gpu")
        assert runtime.assert_gpu_if_requested("cpu") == "cpu"


DEVICE_CHILD = """
import json, sys
import chrombpnet
from jax._src import xla_bridge
from chrombpnet.training import runtime
before = xla_bridge.backends_are_initialized()
result = runtime.assert_gpu_if_requested(sys.argv[1])
state = {"before": before, "result": result, "after": xla_bridge.backends_are_initialized(),
         "clients": sorted(xla_bridge._backends)}
import jax
state["platforms"] = jax.config.jax_platforms
state["default_backend"] = jax.default_backend()
print("STATE " + json.dumps(state))
"""


def device_check_in_fresh_process(device):
    env = {k: v for k, v in os.environ.items() if k != "JAX_PLATFORMS"}
    proc = subprocess.run([sys.executable, "-c", DEVICE_CHILD, device], capture_output=True, text=True, env=env,
                          timeout=300)
    assert proc.returncode == 0, proc.stderr
    line = [x for x in proc.stdout.splitlines() if x.startswith("STATE ")][-1]
    return json.loads(line[len("STATE "):]), proc.stdout


def test_device_cpu_restricts_jax_before_any_backend_exists():
    # --device cpu must not create the CUDA client (a context on every visible GPU): only the CPU one
    state, stdout = device_check_in_fresh_process("cpu")
    assert state["before"] is False
    assert state["result"] == "cpu" and state["platforms"] == "cpu"
    assert state["clients"] == ["cpu"] and state["default_backend"] == "cpu"
    assert "jax backend: cpu" in stdout


def test_device_cpu_after_a_gpu_backend_exists(monkeypatch):
    # too late for jax_platforms: fall back to making the CPU the default device, with a warning
    monkeypatch.setattr(runtime, "backends_initialized", lambda: True)
    monkeypatch.setattr(jax, "default_backend", lambda: "gpu")
    platforms, default_device = jax.config.jax_platforms, jax.config.jax_default_device
    try:
        with pytest.warns(UserWarning, match="initialised already"):
            assert runtime.assert_gpu_if_requested("cpu") == "cpu"
        assert jax.config.jax_default_device == jax.devices("cpu")[0]
        assert jax.config.jax_platforms == platforms
    finally:
        jax.config.update("jax_default_device", default_device)


def test_device_auto_does_not_create_backends():
    # the CLI checks --device before hours of preprocessing: 'auto' leaves the backends to the first model
    state, stdout = device_check_in_fresh_process("auto")
    assert state["before"] is False and state["after"] is False
    assert state["result"] is None and state["platforms"] is None
    assert "first use" in stdout


@pytest.mark.gpu
def test_gpu_present():
    assert runtime.assert_gpu_if_requested("gpu") in ("gpu", "cuda", "rocm")
