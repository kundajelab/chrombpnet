"""Every chrombpnet module imports without TensorFlow and its companions, on the Keras JAX backend.

All imports run in one child interpreter. In that interpreter TensorFlow, tensorflow_probability, tf_keras,
deepdish, deeplift and shap are blocked (`sys.modules[name] = None` makes `import name` raise ImportError), and
KERAS_BACKEND is removed from the environment, so chrombpnet/__init__.py has to select JAX on its own.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PACKAGE_DIR = REPO / "chrombpnet"
SKIP_DIRS = {"figure_notebooks", "notebooks", "__pycache__"}
BLOCKED = ("tensorflow", "tensorflow_probability", "tf_keras", "deepdish", "deeplift", "shap")

# Modules known not to import yet. Each entry becomes a strict xfail, so remove it as soon as the module imports.
# Keep this empty on the main branch.
EXPECTED_IMPORT_FAILURES = set()

# What downstream pipelines that import chrombpnet use from it. None means the module itself.
DOWNSTREAM_CONTRACT = [
    ("chrombpnet.CHROMBPNET", "main"),
    ("chrombpnet.pipelines", None),
    ("chrombpnet.helpers.preprocessing.reads_to_bigwig", "main"),
    ("chrombpnet.evaluation.interpret.interpret", "main"),
    ("chrombpnet.training.predict", "main"),
    ("chrombpnet.training.utils.losses", "multinomial_nll"),
    ("chrombpnet.training.utils.one_hot", "dna_to_one_hot"),
    ("chrombpnet.training.utils.data_utils", "get_cts"),
    ("chrombpnet.helpers.hyperparameters.param_utils", "load_model_wrapper"),
    ("chrombpnet.evaluation.make_bigwigs.bigwig_helper", "get_seq"),
    ("chrombpnet.evaluation.make_bigwigs.bigwig_helper", "read_chrom_sizes"),
    ("chrombpnet.evaluation.make_bigwigs.bigwig_helper", "get_regions"),
    ("chrombpnet.evaluation.make_bigwigs.bigwig_helper", "write_bigwig"),
    ("chrombpnet.helpers.generate_reports.make_html", "main"),
    ("chrombpnet.helpers.generate_reports.make_html_bias", "main"),
    ("chrombpnet.evaluation.modisco.convert_html_to_pdf", "main"),
    ("chrombpnet.data", "DefaultDataFile"),
    ("chrombpnet.data", "get_default_data_path"),
]

CHILD = r"""
import importlib, json, sys, traceback

for name in json.loads(sys.argv[1]):
    sys.modules[name] = None
modules, contract, out = json.loads(sys.argv[2]), json.loads(sys.argv[3]), sys.argv[4]


def describe(exc):
    return "".join(traceback.format_exception_only(type(exc), exc)).strip()


result = {"failures": {}, "contract": {}}
for mod in modules:
    try:
        importlib.import_module(mod)
    except BaseException as e:  # also SystemExit from argparse at import time
        result["failures"][mod] = describe(e)

for mod, attr in contract:
    key = mod if attr is None else mod + ":" + attr
    try:
        obj = importlib.import_module(mod)
        if attr is not None:
            obj = getattr(obj, attr)
            if not callable(obj):
                raise TypeError(key + " is not callable")
        result["contract"][key] = "ok"
    except BaseException as e:
        result["contract"][key] = describe(e)

import keras
result["backend"] = keras.backend.backend()
try:
    from jax._src import xla_bridge
    result["devices_initialised"] = bool(xla_bridge.backends_are_initialized())
except Exception:
    result["devices_initialised"] = None
with open(out, "w") as f:
    json.dump(result, f)
"""


def discover_modules():
    names = []
    for path in sorted(PACKAGE_DIR.rglob("*.py")):
        rel = path.relative_to(REPO)
        if any(part in SKIP_DIRS or part.startswith(".") for part in rel.parts[:-1]):
            continue
        parts = list(rel.with_suffix("").parts)
        if parts[-1] == "__init__":
            parts = parts[:-1]
        names.append(".".join(parts))
    return names


MODULES = discover_modules()


def contract_key(mod, attr):
    return mod if attr is None else mod + ":" + attr


@pytest.fixture(scope="module")
def import_report(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("imports")
    out = tmp / "report.json"
    env = dict(os.environ)
    env.pop("KERAS_BACKEND", None)
    env["KERAS_HOME"] = str(tmp / "keras")
    env["MPLBACKEND"] = "Agg"
    env["PYTHONPATH"] = os.pathsep.join([str(REPO)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    args = [sys.executable, "-c", CHILD, json.dumps(BLOCKED), json.dumps(MODULES),
            json.dumps(DOWNSTREAM_CONTRACT), str(out)]
    proc = subprocess.run(args, cwd=tmp, env=env, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not out.exists():
        pytest.fail("import child exited with {}:\n{}\n{}".format(proc.returncode, proc.stdout[-4000:],
                                                                  proc.stderr[-4000:]))
    with open(out) as f:
        return json.load(f)


def test_discovers_the_package():
    assert "chrombpnet" in MODULES
    assert "chrombpnet.training.utils.model_io" in MODULES
    assert not any("figure_notebooks" in m for m in MODULES)


@pytest.mark.parametrize("module", [
    pytest.param(m, marks=pytest.mark.xfail(strict=True, reason="listed in EXPECTED_IMPORT_FAILURES"))
    if m in EXPECTED_IMPORT_FAILURES else m
    for m in MODULES
])
def test_module_imports_without_tensorflow(import_report, module):
    error = import_report["failures"].get(module)
    assert error is None, "{} failed to import: {}".format(module, error)


def test_keras_backend_is_jax(import_report):
    assert import_report["backend"] == "jax"


@pytest.mark.parametrize("mod,attr", DOWNSTREAM_CONTRACT, ids=[contract_key(m, a) for m, a in DOWNSTREAM_CONTRACT])
def test_downstream_contract(import_report, mod, attr):
    status = import_report["contract"][contract_key(mod, attr)]
    assert status == "ok", status


def test_imports_do_not_initialise_jax_devices(import_report):
    # Importing chrombpnet must not create a JAX client (on a shared GPU that grabs a CUDA context at import time).
    if import_report["devices_initialised"] is None:
        pytest.skip("this jax version does not expose backends_are_initialized")
    assert import_report["devices_initialised"] is False
