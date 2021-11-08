import tensorflow as tf
import tensorflow_probability as tfp
import sys
import pdb 
import numpy as np 


#from https://github.com/kundajelab/basepair/blob/cda0875571066343cdf90aed031f7c51714d991a/basepair/losses.py#L87
def multinomial_nll(true_counts, logits):
    """Compute the multinomial negative log-likelihood
    Args:
      true_counts: observed count values
      logits: predicted logit values
    """
    counts_per_example = tf.reduce_sum(true_counts, axis=-1)
    dist = tfp.distributions.Multinomial(total_count=counts_per_example,
                                         logits=logits)
    return (-tf.reduce_sum(dist.log_prob(true_counts)) / 
            tf.cast(tf.shape(true_counts)[0], dtype=tf.float32))


class MultichannelMultinomialNLL(object):
    def __init__(self, n, weights=None):
        self.__name__ = "MultichannelMultinomialNLL"
        self.n = n
        if weights is None:
            self.weights = [1]*self.n
        else:
            self.weights = weights
        #print(self.weights, self.n)


    def __call__(self, true_counts, logits):
        for i in range(self.n):
            print(true_counts[..., i], logits[..., i])
            loss = multinomial_nll(true_counts[..., i], logits[..., i])
            #print(loss)
            if i == 0:
                total = self.weights[i]*loss
            else:
                total += self.weights[i]*loss
        return total

    def get_config(self):
        return {"n": self.n, "weights":self.weights}


