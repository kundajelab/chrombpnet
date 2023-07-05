import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.backend import int_shape

def bpnet_seq(inputlen=2114, outputlen=1000, filters=512, ndil=8):
    """
    Classic BPNet architecture with sequence-only input. Predicts profile
    logits and log-counts.

    Inputs and outputs have length inputlen and outputlen respectively.
    """

    inp = keras.Input((inputlen,4))

    x = keras.layers.Conv1D(filters,
                            kernel_size=21,
                            padding='valid',
                            activation='relu')(inp)

    for i in range(1, ndil+1):
        conv_x = keras.layers.Conv1D(filters,
                                     kernel_size=3,
                                     padding='valid',
                                     activation='relu',
                                     dilation_rate=2**i)(x)

        x_len = int_shape(x)[1]
        conv_x_len = int_shape(conv_x)[1]

        assert((x_len - conv_x_len) % 2 == 0) # Necessary for symmetric cropping

        x = keras.layers.Cropping1D((x_len - conv_x_len) // 2)(x)
        x = keras.layers.add([conv_x, x])

    prof = keras.layers.Conv1D(1, 75, padding='valid')(x)
    prof = keras.layers.Flatten(name="logits")(prof)
    assert prof.shape[1] == outputlen, prof.shape[1]

    ct = keras.layers.GlobalAvgPool1D()(x)
    ct = keras.layers.Dense(1, name="logcounts")(ct)

    return keras.Model(inputs=inp, outputs=[prof,ct])


def chrombpnet(inputlen=2114, outputlen=1000, filters=512, ndil=8):
    """
    ChromBPNet architecture with sequence + bias counts + bias profile logits
    as inputs. Passes sequence through a bpnet_seq module. The counts output
    are logsumexp-ed with bias counts input, and the profile output is simply
    added to the bias logits input.
    """

    # Inputs: sequence + (predicted) bias logits + (predicted) bias logcounts
    inp = keras.Input((inputlen,4))
    inp_bias_logits = keras.Input(shape=(outputlen,))
    inp_bias_logcounts = keras.Input(shape=(1,))

    # sequence is processed by bpnet_seq model
    seq_model = bpnet_seq(inputlen=inputlen,
                          outputlen=outputlen,
                          filters=filters,
                          ndil=ndil)
    unbias_logits, unbias_logcounts = seq_model(inp)

    # PROFILE: predicted logits are simply a sum of unbias logits with
    # (predicted) bias logits
    out_logits = keras.layers.Add(name="logits_w_bias")([unbias_logits, inp_bias_logits]) # <== add in bias logits

    # COUNTS: final count is log(exp(unbias_logcounts) + exp(inp_bias_logcounts)))
    concat_cts = keras.layers.concatenate([unbias_logcounts, inp_bias_logcounts],
                                          axis=-1)
    out_logcounts = keras.layers.Lambda(
                        lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True),
                        name="logcounts_w_bias")(concat_cts)

    return keras.Model(inputs=[inp, inp_bias_logits, inp_bias_logcounts],
                       outputs=[out_logits, out_logcounts])
