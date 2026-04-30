import numpy as np
import re
from tensorflow.keras.layers import Input, Conv1D, Dense, Concatenate, Flatten, Reshape
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.utils import get_custom_objects
from chrombpnet.training.utils.losses import multinomial_nll
import tensorflow as tf
import random as rn
import os

os.environ['PYTHONHASHSEED'] = '0'

# Profile precrop conv and counts Dense produced by bpnet_model.py
# (and the wo_bias_* variants used by chrombpnet_with_bias_model.py's trunk).
_TRAINABLE_HEAD_NAMES = {
    "logcount_predictions", "logcounts",
    "prof_out_precrop", "wo_bias_bpnet_prof_out_precrop",
}

# Matches the dilated conv layers in bpnet_model.py: `bpnet_{i}conv` for
# i=1..n_dil_layers (and the wo_bias_ variant). The deepest one (highest
# index) is the last dilated conv, which we also keep trainable.
_DILATED_CONV_RE = re.compile(r"^(?:wo_bias_)?bpnet_(\d+)conv$")


def _load_pretrained_bias(path):
    get_custom_objects().update({"multinomial_nll": multinomial_nll, "tf": tf})
    return load_model(path, compile=False)


def _last_dilated_conv_name(submodel):
    indexed = []
    for L in submodel.layers:
        m = _DILATED_CONV_RE.match(L.name)
        if m:
            indexed.append((int(m.group(1)), L.name))
    assert indexed, (
        "no `bpnet_{i}conv` (or wo_bias_ variant) layers found in bias "
        "submodel — cannot identify the last dilated conv to fine-tune"
    )
    indexed.sort()
    return indexed[-1][1]


def _freeze_except_heads_and_last_dilated(submodel):
    last_dil = _last_dilated_conv_name(submodel)
    trainable_set = _TRAINABLE_HEAD_NAMES | {last_dil}
    n_unfrozen = 0
    for L in submodel.layers:
        if L.name in trainable_set:
            L.trainable = True
            n_unfrozen += 1
        else:
            L.trainable = False
    return n_unfrozen


def getModelGivenModelOptionsAndWeightInits(args, model_params):
    counts_loss_weight = float(model_params['counts_loss_weight'])
    sequence_len = int(model_params['inputlen'])
    out_pred_len = int(model_params['outputlen'])
    paths = [p for p in model_params['pretrained_bias_model_paths'].split(',') if p]
    assert len(paths) >= 1, "pretrained_bias_model_paths is empty"
    combiner_kernel = int(model_params.get('combiner_profile_kernel_size', 1))

    seed = args.seed
    np.random.seed(seed)
    tf.random.set_seed(seed)
    rn.seed(seed)

    inp = Input(shape=(sequence_len, 4), name='sequence')

    profile_branches = []
    counts_branches = []
    for i, p in enumerate(paths):
        sub = _load_pretrained_bias(p)
        sub._name = "pretrained_bias_{}".format(i)
        n_unfrozen = _freeze_except_heads_and_last_dilated(sub)
        assert n_unfrozen == 3, (
            "pretrained bias model {} has {} unfrozen layers, expected 3 "
            "(last dilated conv `bpnet_{{n}}conv`, profile precrop "
            "`prof_out_precrop`, counts dense `logcount_predictions` — or "
            "their wo_bias_ variants). Path: {}".format(i, n_unfrozen, p)
        )
        prof_i, lc_i = sub(inp)
        prof_i = Reshape((out_pred_len, 1), name="bias{}_profile_reshape".format(i))(prof_i)
        profile_branches.append(prof_i)
        counts_branches.append(lc_i)

    if len(profile_branches) > 1:
        profile_concat = Concatenate(axis=-1, name="profile_combiner_concat")(profile_branches)
        counts_concat = Concatenate(axis=-1, name="counts_combiner_concat")(counts_branches)
    else:
        profile_concat = profile_branches[0]
        counts_concat = counts_branches[0]

    profile_combined = Conv1D(
        filters=1, kernel_size=combiner_kernel, padding='same',
        name="profile_combiner_conv",
    )(profile_concat)
    profile_out = Flatten(name="logits_profile_predictions")(profile_combined)

    # Dense is created last so the parent model has a clearly identifiable
    # counts head. Named `combined_logcount_predictions` (not the bare
    # `logcount_predictions`) so its weights don't collide with the same-named
    # Dense layer inside each loaded sub-model when Keras writes the H5 — that
    # collision otherwise raises "Unable to synchronously create dataset
    # (name already exists)" at the first ModelCheckpoint save.
    # `find_chrombpnet_hyperparams.adjust_bias_model_logcounts` knows to look
    # up this name (in addition to the original two) so the post-finetune
    # δ-scale step still finds the right Dense.
    count_out = Dense(1, name="combined_logcount_predictions")(counts_concat)

    model = Model(inputs=[inp], outputs=[profile_out, count_out])

    # 3 trainable layers per bias model (last dilated conv + profile precrop
    # conv + counts Dense), each contributing kernel + bias = 2 weight tensors.
    # Plus combiner Conv1D (2 tensors) and Dense (2 tensors).
    expected_n_trainable = 2 * 3 * len(paths) + 2 + 2
    got = len(model.trainable_weights)
    assert got == expected_n_trainable, (
        "expected {} trainable weight tensors, got {}".format(expected_n_trainable, got)
    )

    model.compile(
        optimizer=Adam(learning_rate=args.learning_rate),
        loss=[multinomial_nll, 'mse'],
        loss_weights=[1, counts_loss_weight],
    )
    return model


def save_model_without_bias(model, output_prefix):
    # The fine-tuned combiner *is* the bias model. The trained weights are
    # already saved by train.fit_and_evaluate via model.save(...). Nothing
    # additional to extract here.
    return
