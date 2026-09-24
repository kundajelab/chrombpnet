import os

__version__ = "2.0.0.dev0"

# ChromBPNet runs on Keras 3 with the JAX backend (TensorFlow has no CUDA 13 build). This must happen before
# anything imports keras; an explicit KERAS_BACKEND from the user still wins.
os.environ.setdefault("KERAS_BACKEND", "jax")
# Allocate GPU memory on demand instead of reserving 75% of the device at start-up (the pipeline keeps the GPU
# for hours and often shares it with other jobs). Set XLA_PYTHON_CLIENT_MEM_FRACTION to cap it.
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
