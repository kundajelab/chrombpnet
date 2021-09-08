import tensorflow as tf
import tensorflow_probability as tfp
import sys
import pdb 
import numpy as np 
def get_weighted_binary_crossentropy(w0_weights, w1_weights,ambig_val=np.nan):
    # Compute the task-weighted cross-entropy loss, where every task is weighted by 1 - (fraction of non-ambiguous examples that are positive)
    # In addition, weight everything with label -1 to 0
    import numpy as np 
    w0_weights=np.array(w0_weights);
    w1_weights=np.array(w1_weights);
    def weighted_binary_crossentropy(y_true,y_pred):
        weightsPerTaskRep = y_true*w1_weights[None,:] + (1-y_true)*w0_weights[None,:]
        nonAmbig=tf.math.logical_not(tf.math.is_nan(y_true))
        nonAmbigTimesWeightsPerTask = tf.boolean_mask(weightsPerTaskRep,nonAmbig)
        return tf.mean(tf.binary_crossentropy(tf.boolean_mask(y_true,nonAmbig),tf.boolean_mask(y_pred,nonAmbig))*nonAmbigTimesWeightsPerTask, axis=-1);
    return weighted_binary_crossentropy; 

def ambig_binary_crossentropy(y_true,y_pred):
    nonAmbig=tf.math.logical_not(tf.math.is_nan(y_true))
    return tf.mean(tf.binary_crossentropy(tf.boolean_mask(y_true,nonAmbig), tf.boolean_mask(y_pred,nonAmbig)), axis=-1);

def ambig_mean_squared_error(y_true, y_pred):
    nonAmbig=tf.math.logical_not(tf.math.is_nan(y_true))
    return tf.mean(tf.square(tf.boolean_mask(y_pred,nonAmbig) - tf.boolean_mask(y_true,nonAmbig)), axis=-1)

def ambig_mean_absolute_error(y_true, y_pred):
    nonAmbig=tf.math.logical_not(tf.math.is_nan(y_true))
    return tf.mean(tf.abs(tf.boolean_mask(y_pred,nonAmbig) - tf.boolean_mask(y_true,nonAmbig)), axis=-1)

def ambig_log_poisson(y_true,y_pred):
    nonAmbig=tf.math.logical_not(tf.math.is_nan(y_true))
    return tf.nn.log_poisson_loss(nonAmbig,tf.log(y_pred),compute_full_loss=True)


#PROFILE MODEL LOSSES #
def get_loss_weights(tdb_path,chrom,label_attribute,ambig_attribute,upsample_attribute,tdb_partition_thresh_for_upsample):
    import tiledb
    from kerasAC.tiledb_config import get_default_config
    import pdb 
    tdb_config=get_default_config()
    ctx=tiledb.Ctx(tdb_config)
    tdb_array=tiledb.DenseArray(tdb_path+"."+chrom,mode='r',ctx=ctx)
    print("opened:"+tdb_path+"."+chrom+" for reading")
    vals=tdb_array[:]
    print("got tdb vals")
    #label_attribute

def multinomial_nll_probs(true_counts, logits):
    """Compute the multinomial negative log-likelihood
    Args:
      true_counts: observed count values
      logits: predicted logit values
    """
    counts_per_example = tf.reduce_sum(true_counts, axis=-1)
    dist = tfp.distributions.Multinomial(total_count=counts_per_example,
                                         probs=tf.linalg.normalize(logits+1e-10, ord=1, axis=1)[0])
    return (-tf.reduce_sum(dist.log_prob(true_counts)) / 
            tf.cast(tf.shape(true_counts)[0], dtype=tf.float32))



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



def multinomial_nll_masked(true_counts, logits):
    """Compute the multinomial negative log-likelihood
    Args:
      true_counts: observed count values
      logits: predicted logit values
    """
    counts_per_example = tf.reduce_sum(true_counts, axis=-1)
    mask = counts_per_example >= 4.6

    counts_per_example1 = tf.boolean_mask(counts_per_example,mask)
    logits1 = tf.boolean_mask(logits,mask,axis=0)
    true_counts1 = tf.boolean_mask(true_counts,mask,axis=0)

    dist = tfp.distributions.Multinomial(total_count=counts_per_example1,
                                         logits=logits1)
    #print(mask)
    return (-tf.reduce_sum(dist.log_prob(true_counts1)) / 
            tf.cast(tf.shape(true_counts1)[0], dtype=tf.float32))

def custom_mse(y_true, y_pred):
    # calculating squared difference between target and predicted values 
    loss = tf.square(y_pred - y_true)  # (batch_size, 2)
    # summing both loss values along batch dimension 
    loss = tf.sum(loss, axis=1)        # (batch_size,)
    return loss

#from https://github.com/kundajelab/basepair/blob/cda0875571066343cdf90aed031f7c51714d991a/basepair/losses.py#L87
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
            loss = multinomial_nll(true_counts[..., i], logits[..., i])
            #print(loss)
            if i == 0:
                total = self.weights[i]*loss
            else:
                total += self.weights[i]*loss
        return total

    def get_config(self):
        return {"n": self.n, "weights":self.weights}


class MultichannelMultinomialNLLProbs(object):
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
            loss = multinomial_nll_probs(true_counts[..., i], logits[..., i])
            #print(loss)
            if i == 0:
                total = self.weights[i]*loss
            else:
                total += self.weights[i]*loss
        return total

    def get_config(self):
        return {"n": self.n, "weights":self.weights}



class MultichannelMultinomialMaskedNLL(object):
    def __init__(self, n, weights=None):
        self.__name__ = "MultichannelMultinomialMaskedNLL"
        self.n = n
        if weights is None:
            self.weights = [1]*self.n
        else:
            self.weights = weights
        #print(self.weights, self.n)


    def __call__(self, true_counts, logits):
        for i in range(self.n):
            loss = multinomial_nll_masked(true_counts[..., i], logits[..., i])
            #print(loss)
            if i == 0:
                total = self.weights[i]*loss
            else:
                total += self.weights[i]*loss
        return total

    def get_config(self):
        return {"n": self.n, "weights":self.weights}

class MultichannelMultinomialMSE(object):
    def __init__(self, n, weights=None):
        self.__name__ = "MultichannelMultinomialMSE"
        self.n = n
        if weights is None:
            self.weights = [1]*self.n
        else:
            self.weights = weights
        #print(self.weights, self.n)
        self.mse = tf.keras.losses.MeanSquaredError()

    def __call__(self, true_counts, logits):
        for i in range(self.n):
            #print(true_counts[..., i], logits[..., i])
            loss = self.mse(true_counts[..., i], logits[..., i])
            #print(loss, self.n)
            if i == 0:
                total = self.weights[i]*loss
            else:
                total += self.weights[i]*loss
        return total

    def get_config(self):
        return {"n": self.n, "weights":self.weights}


class MultiChannelKLDivergence(object):
    def __init__(self, n, weights=None):
        self.__name__ = "MultiChannelKLDivergence"
        self.n = n
        if weights is None:
            self.weights = [1]*self.n
        else:
            self.weights = weights
        #print(self.weights, self.n)
        #self.cce = tf.keras.losses.CategoricalCrossentropy(from_logits=True)
        self.cce = tf.keras.losses.KLDivergence()

    def __call__(self, true_counts, logits):
        for i in range(self.n):
            #tf.print("true_counts",true_counts, output_stream=sys.stdout)
            #tf.print("logits",logits,output_stream=sys.stdout)

            loss = self.cce(true_counts[..., i], logits[..., i])

            #print(loss, self.n)
            if i == 0:
                total = self.weights[i]*loss
            else:
                total += self.weights[i]*loss
        return total

    def get_config(self):
        return {"n": self.n, "weights":self.weights}

