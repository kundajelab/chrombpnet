"""Optimizer factory for the architecture files (--optimizer, --muon-lr, --ema, --lr-schedule).

The defaults reproduce chrombpnet 1.x: Adam(learning_rate=args.learning_rate), constant learning rate, no EMA.
"""
import re

import keras
from keras import ops

DEFAULT_MUON_LR = 2e-3
# The dilated residual convs (bpnet_1conv ... bpnet_Nconv, with or without the wo_bias_ prefix). The first conv
# (bpnet_1st_conv), the profile head (prof_out_precrop), the count Dense and all biases stay on Adam, following
# the reference Muon recipe (hidden-to-hidden weights only).
MUON_VARIABLES = r"bpnet_\d+conv/kernel"


@keras.saving.register_keras_serializable(package="chrombpnet")
class ChromBPNetMuon(keras.optimizers.Muon):
    """keras.optimizers.Muon for dilated Conv1D kernels.

    Stock Muon sends every non-2D variable to its AdamW path, so on BPNet (all kernels are 3D Conv1D kernels) it
    is just AdamW. Here the kernels whose path matches `muon_variables` get the Muon update: the (k, in, out)
    momentum is flattened to (k*in, out) for Newton-Schulz and the rms_rate scaling, then reshaped back. Every
    other variable uses the Adam path at learning_rate * adam_lr_ratio.
    """

    def __init__(self, muon_variables=MUON_VARIABLES, **kwargs):
        super().__init__(**kwargs)
        self.muon_variables = muon_variables

    def _should_use_adamw(self, variable):
        if len(variable.shape) not in (2, 3) or re.search(self.muon_variables, variable.path) is None:
            return True
        return any(re.search(keyword, variable.path) for keyword in self.exclude_layers)

    def _muon_update_step(self, gradient, variable, lr, m):
        self.assign_add(m, ops.add(gradient, m * (self.momentum - 1)))
        if self.nesterov:
            g = ops.add(gradient, self.momentum * m)
        else:
            g = m
        shape = variable.shape
        g = ops.reshape(g, (-1, shape[-1]))
        update = self.lr_adjust(lr * self.zeropower_via_newtonschulz5(g, self.ns_steps))
        self.assign_sub(variable, ops.reshape(update, shape))

    def get_config(self):
        config = super().get_config()
        config["muon_variables"] = self.muon_variables
        return config


def _arg(args, name, default):
    value = getattr(args, name, None)
    return default if value is None else value


def make_learning_rate(peak, args):
    """A float (constant) or a warmup + cosine schedule.

    --lr-schedule cosine reads args.total_steps (required; train.py fills it with epochs * batches per epoch) and
    args.warmup_steps (default 10% of total_steps): linear warmup from 1% of the peak to the peak, then cosine
    decay to args.lr_min_ratio (default 0.01) of the peak at total_steps. The schedule assumes the run lasts
    total_steps; early stopping (--early-stop) may end it during the decay.
    """
    schedule = _arg(args, "lr_schedule", "constant")
    if schedule == "constant":
        return peak
    if schedule != "cosine":
        raise ValueError("lr_schedule must be 'constant' or 'cosine', got {!r}".format(schedule))
    total_steps = _arg(args, "total_steps", 0)
    if not total_steps:
        raise ValueError("--lr-schedule cosine needs args.total_steps (the number of training steps)")
    total_steps = int(total_steps)
    warmup_steps = int(_arg(args, "warmup_steps", total_steps // 10))
    if not 0 <= warmup_steps < total_steps:
        raise ValueError("warmup_steps must be in [0, total_steps), got {} for total_steps={}".format(
            warmup_steps, total_steps))
    alpha = float(_arg(args, "lr_min_ratio", 0.01))
    if warmup_steps == 0:
        return keras.optimizers.schedules.CosineDecay(peak, decay_steps=total_steps, alpha=alpha)
    return keras.optimizers.schedules.CosineDecay(0.01 * peak, decay_steps=total_steps - warmup_steps, alpha=alpha,
                                                  warmup_target=peak, warmup_steps=warmup_steps)


def make_optimizer(args):
    """Optimizer from args (all attributes optional, read with getattr):

    optimizer      'adam' (default) or 'muon'
    learning_rate  Adam learning rate (default 1e-3); with muon, the learning rate of the Adam-routed variables
    muon_lr        learning rate of the Muon-routed conv kernels (default 2e-3)
    muon_weight_decay  decoupled weight decay of the Muon-routed kernels (default 0: none)
    ema            keep an exponential moving average of the weights (default False); train.py then validates
                   and saves the EMA weights (SwapEMAWeights)
    ema_momentum   default 0.999
    lr_schedule    'constant' (default) or 'cosine' (see make_learning_rate)
    """
    name = _arg(args, "optimizer", "adam").lower()
    learning_rate = float(_arg(args, "learning_rate", 1e-3))
    ema = dict(use_ema=bool(_arg(args, "ema", False)), ema_momentum=float(_arg(args, "ema_momentum", 0.999)))
    if name == "adam":
        return keras.optimizers.Adam(learning_rate=make_learning_rate(learning_rate, args), **ema)
    if name == "muon":
        muon_lr = float(_arg(args, "muon_lr", DEFAULT_MUON_LR))
        weight_decay = float(_arg(args, "muon_weight_decay", 0.0))
        return ChromBPNetMuon(learning_rate=make_learning_rate(muon_lr, args), adam_lr_ratio=learning_rate / muon_lr,
                              adam_weight_decay=None, weight_decay=weight_decay or None, exclude_embeddings=False,
                              momentum=0.95, nesterov=True, ns_steps=5, rms_rate=0.2, **ema)
    raise ValueError("optimizer must be 'adam' or 'muon', got {!r}".format(name))
