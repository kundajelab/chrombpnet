# syntax=docker/dockerfile:1
# ChromBPNet 2.x (Keras 3 on JAX). CUDA and cuDNN come from JAX's pip wheels inside the pixi environment, so the
# image needs no CUDA base: the host only provides the NVIDIA driver (>= 580 for cuda13, H100 / B200 / RTX PRO
# 6000 Blackwell), injected by the NVIDIA Container Toolkit (`docker run --gpus all`) or `apptainer --nv`.
#
#   docker build -t chrombpnet:cuda13 .
#   docker build --build-arg PIXI_ENV=cuda12 --build-arg CUDA_VERSION=12.0 -t chrombpnet:cuda12 .
#   docker build --build-arg PIXI_ENV=default -t chrombpnet:cpu .
#   docker run --rm --gpus all chrombpnet:cuda13 python -c "import jax; print(jax.devices())"
#   apptainer run --nv --cleanenv docker://ghcr.io/kundajelab/chrombpnet:<tag> chrombpnet --help
ARG PIXI_VERSION=0.81.0

FROM ghcr.io/prefix-dev/pixi:${PIXI_VERSION}-trixie-slim AS build
ARG PIXI_ENV=cuda13
# The cuda envs require the __cuda virtual package and the builder has no GPU driver. Build time only.
ARG CUDA_VERSION=13.0
WORKDIR /opt/chrombpnet
# Third-party packages first: this layer is reused until the manifest or the lock file changes.
COPY pyproject.toml pixi.lock conda-pypi-map.json ./
RUN CONDA_OVERRIDE_CUDA=${CUDA_VERSION} pixi install --locked -e ${PIXI_ENV} --skip chrombpnet
COPY README.md LICENSE MANIFEST.in ./
COPY chrombpnet ./chrombpnet
RUN CONDA_OVERRIDE_CUDA=${CUDA_VERSION} pixi install --locked -e ${PIXI_ENV} \
 && pixi shell-hook -e ${PIXI_ENV} --shell=bash --as-is > /opt/chrombpnet/shell-hook.sh \
 && printf '#!/bin/bash\n. /opt/chrombpnet/shell-hook.sh\nexec "$@"\n' > /opt/chrombpnet/entrypoint.sh \
 && chmod 0755 /opt/chrombpnet/entrypoint.sh \
 && CONDA_OVERRIDE_CUDA=${CUDA_VERSION} JAX_PLATFORMS=cpu pixi run -e ${PIXI_ENV} --as-is python -c \
    "import chrombpnet, keras, jax; print('chrombpnet', chrombpnet.__version__, '| keras', keras.__version__, keras.backend.backend(), '| jax', jax.__version__)"

FROM debian:trixie-slim AS runtime
ARG PIXI_ENV=cuda13
# Optional Google Cloud CLI (gsutil/gcloud) for pipelines that stage data in buckets; chrombpnet does not need it.
ARG WITH_GCLOUD=false
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates procps \
 && if [ "${WITH_GCLOUD}" = "true" ]; then \
      apt-get install -y --no-install-recommends curl gnupg \
      && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
      && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
         > /etc/apt/sources.list.d/google-cloud-sdk.list \
      && apt-get update && apt-get install -y --no-install-recommends google-cloud-cli; \
    fi \
 && rm -rf /var/lib/apt/lists/*
# Same prefix as in the build stage: conda environments are not relocatable.
COPY --from=build /opt/chrombpnet /opt/chrombpnet
# The entrypoint sources the pixi activation; these also cover `apptainer exec` and `docker run --entrypoint`,
# which bypass it. No LD_LIBRARY_PATH: it would shadow the CUDA libraries of the JAX wheels.
# Apptainer passes the host environment and binds $HOME: ignore the user site (~/.local/lib/python3.x) and any
# host PYTHONPATH/PYTHONHOME (e.g. from `module load`), so only the image's packages are imported. Image values
# win over host ones (Apptainer unsets a host variable that the image sets to empty); Python treats empty as unset.
ENV PATH=/opt/chrombpnet/.pixi/envs/${PIXI_ENV}/bin:${PATH} \
    CONDA_PREFIX=/opt/chrombpnet/.pixi/envs/${PIXI_ENV} \
    PYTHONNOUSERSITE=1 \
    PYTHONPATH="" \
    PYTHONHOME="" \
    KERAS_BACKEND=jax \
    XLA_PYTHON_CLIENT_PREALLOCATE=false \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    KERAS_HOME=/tmp/keras \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8
WORKDIR /work
ENTRYPOINT ["/opt/chrombpnet/entrypoint.sh"]
CMD ["chrombpnet", "--help"]
