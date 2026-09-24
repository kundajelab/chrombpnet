"""Static checks of the infrastructure files: pixi activation, CI action refs and the Docker image environment.

These need no network and no conda tools, so they also run in the uv-only CI job.
"""
import re
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
WORKFLOWS = REPO / ".github" / "workflows"
DOCKERFILE = REPO / "Dockerfile"

# Actions that publish only full version tags (setup-uv has no floating `v10` tag or branch since v8): a
# major-only ref fails at "Set up job" before any step runs.
FULL_VERSION_ONLY = {"astral-sh/setup-uv"}
FULL_VERSION = re.compile(r"v\d+\.\d+\.\d+")
COMMIT_SHA = re.compile(r"[0-9a-f]{40}")


def _require(path):
    if not path.exists():
        pytest.skip("{} is not in this checkout".format(path.relative_to(REPO)))
    return path


def _uses(path):
    return re.findall(r"^\s*-?\s*uses:\s*([^\s#]+)", path.read_text(), flags=re.MULTILINE)


def _runtime_env(dockerfile_text):
    """KEY=VALUE pairs of the ENV instructions of the last stage (the runtime image)."""
    stage = dockerfile_text.split("\nFROM ")[-1]
    env = {}
    for block in re.findall(r"^ENV (.*?)(?<!\\)$", stage, flags=re.MULTILINE | re.DOTALL):
        for key, value in re.findall(r'([A-Za-z_][A-Za-z0-9_]*)=("[^"]*"|\S*)', block.replace("\\\n", " ")):
            env[key] = value.strip('"')
    return env


def test_pixi_activation_leaves_preallocation_to_the_caller():
    # pixi activation variables override what the caller exported, so the on-demand allocation default lives in
    # chrombpnet/__init__.py (setdefault) and every task that touches the GPU imports chrombpnet first.
    config = tomllib.loads(_require(REPO / "pyproject.toml").read_text())
    pixi = config["tool"]["pixi"]
    assert "XLA_PYTHON_CLIENT_PREALLOCATE" not in pixi["activation"]["env"]
    assert pixi["activation"]["env"]["KERAS_BACKEND"] == "jax"
    for feature in ("cuda13", "cuda12"):
        assert "import chrombpnet" in pixi["feature"][feature]["tasks"]["gpu-check"]


def test_workflow_action_refs_are_pinned():
    workflows = sorted(_require(WORKFLOWS).glob("*.yml"))
    assert workflows
    for path in workflows:
        for use in _uses(path):
            action, sep, ref = use.partition("@")
            assert sep and ref, "{}: `uses: {}` has no ref".format(path.name, use)
            if action in FULL_VERSION_ONLY:
                assert FULL_VERSION.fullmatch(ref) or COMMIT_SHA.fullmatch(ref), (
                    "{}: {} publishes no major-only tags; pin a full version such as v1.2.3, not @{}".format(
                        path.name, action, ref))


def test_workflow_uv_jobs_use_setup_uv():
    uses = _uses(_require(WORKFLOWS / "ci.yml"))
    assert [u for u in uses if u.startswith("astral-sh/setup-uv@")], "ci.yml no longer installs uv"


def test_runtime_image_isolates_python_from_host():
    # Apptainer passes the host environment and binds $HOME; image ENV values win over host ones.
    env = _runtime_env(_require(DOCKERFILE).read_text())
    assert env.get("PYTHONNOUSERSITE") == "1"
    assert env.get("PYTHONPATH") == ""
    assert env.get("PYTHONHOME") == ""
    assert env.get("KERAS_BACKEND") == "jax"
    assert env.get("XLA_PYTHON_CLIENT_PREALLOCATE") == "false"


def test_runtime_env_parser():
    text = 'FROM a AS build\nENV X=1\nFROM b\nENV A=1 \\\n    B="" \\\n    C=x,y\nENV D=\n'
    assert _runtime_env(text) == {"A": "1", "B": "", "C": "x,y", "D": ""}
