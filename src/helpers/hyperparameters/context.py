import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../training/')))
import utils.losses as losses
import utils.one_hot as one_hot

import tensorflow as tf
import tensorflow_probability as tfp

def multinomial_nll(true_counts, logits):
    """Compute the multinomial negative log-likelihood
    Args:
      true_counts: observed count values [B x L]
      logits: predicted logit values [B x L]
    """
    counts_per_example = tf.reduce_sum(true_counts, axis=-1)
    dist = tfp.distributions.Multinomial(total_count=counts_per_example,
                                         logits=logits)
    return (-tf.reduce_sum(dist.log_prob(true_counts)) / 
            tf.cast(tf.shape(true_counts)[0], dtype=tf.float32))

def load_model_wrapper(model_h5):
    # read .h5 model
    import tensorflow as tf
    from keras.utils.generic_utils import get_custom_objects
    from tensorflow.keras.models import load_model
    custom_objects={"MultichannelMultinomialNLL": losses.MultichannelMultinomialNLL, "tf": tf, "multinomial_nll":multinomial_nll}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_h5)
    print("got the model")
    model.summary()
    return model

