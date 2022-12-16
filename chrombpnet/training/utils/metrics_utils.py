
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from scipy.stats import multinomial
import math
from scipy.spatial.distance import jensenshannon


plt.rcParams["figure.figsize"]=10,5
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}
matplotlib.rc('font', **font)


def _fix_sum_to_one(probs):
    """
      Fix probability arrays whose sum is fractinally above or 
      below 1.0
      
      Args:
          probs (numpy.ndarray): An array whose sum is almost equal
              to 1.0
              
      Returns:
          np.ndarray: array that sums to 1
    """
    
    _probs = np.copy(probs)
    
    if np.sum(_probs) > 1.0:        
        _probs[np.argmax(_probs)] -= np.sum(_probs) - 1.0    
    
    if np.sum(_probs) < 1.0:
        _probs[np.argmin(_probs)] += 1.0 - np.sum(_probs)

    return _probs
    

def density_scatter(x, y, xlab, ylab, ax = None, sort = True, bins = 20):
    """
    Scatter plot colored by 2d histogram
    """
    bad_indices=np.where(np.isnan(x))+np.where(np.isnan(y))
    x=x[~np.isin(np.arange(x.size),bad_indices)]
    y=y[~np.isin(np.arange(y.size),bad_indices)]

    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0
    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Density')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    return ax
 
#https://github.com/kundajelab/basepairmodels/blob/cf8e346e9df1bad9e55bd459041976b41207e6e5/basepairmodels/cli/metrics.py#L18
# replacing TracebackExceptions with assertions
def mnll(true_counts, logits=None, probs=None):
    """
        Compute the multinomial negative log-likelihood between true
        counts and predicted values of a BPNet-like profile model
        
        One of `logits` or `probs` must be given. If both are
        given `logits` takes preference.
        Args:
            true_counts (numpy.array): observed counts values
            
            logits (numpy.array): predicted logits values
            
            probs (numpy.array): predicted values as probabilities
          
        Returns:
            float: cross entropy
    
    """

    dist = None 
    
    if logits is not None:
        
        # check for length mismatch
        assert (len(logits) == len(true_counts))
    
        # convert logits to softmax probabilities
        probs = logits - logsumexp(logits)
        probs = np.exp(probs)
        
    elif probs is not None:      
        
        # check for length mistmatch
        assert(len(probs) == len(true_counts))
    
        # check if probs sums to 1
        # why is this nans sometimes
        assert( abs(1.0 - np.sum(probs)) < 1e-1)
         
    else:
        
        # both 'probs' and 'logits' are None
         print(
            "At least one of probs or logits must be provided. "
            "Both are None.")
  
    # compute the nmultinomial distribution
    mnom = multinomial(np.sum(true_counts), probs)
    return -(mnom.logpmf(true_counts) / len(true_counts))
    
#https://github.com/kundajelab/basepairmodels/blob/cf8e346e9df1bad9e55bd459041976b41207e6e5/basepairmodels/cli/metrics.py#L129
def get_min_max_normalized_value(val, minimum, maximum):
    ret_val = (val - maximum) / (minimum - maximum)

    if ret_val < 0:
        return 0
    
    if ret_val > 1:
        return 1
    return ret_val

#https://github.com/kundajelab/basepairmodels/blob/cf8e346e9df1bad9e55bd459041976b41207e6e5/basepairmodels/cli/fastpredict.py#L59
def mnll_min_max_bounds(profile):
    """
        Min Max bounds for the mnll metric
        
        Args:
            profile (numpy.ndarray): true profile 
        Returns:
            tuple: (min, max) bounds values
    """
    
    # uniform distribution profile
    uniform_profile = np.ones(len(profile)) * (1.0 / len(profile))

    # profile as probabilities
    profile = profile.astype(np.float64)
    
    # profile as probabilities
    profile_prob = profile / np.sum(profile)
    
    # the scipy.stats.multinomial function is very sensitive to 
    # profile_prob summing to exactly 1.0, if not you get NaN as the
    # resuls. In majority of the cases we can fix that problem by
    # adding or substracting the difference (but unfortunately it
    # doesnt always and there are cases where we still see NaNs, and
    # those we'll set to 0)
    profile_prob = _fix_sum_to_one(profile_prob)
    #print(profile, profile_prob)

    # mnll of profile with itself
    min_mnll = mnll(profile, probs=profile_prob)
    
    # if we still find a NaN, even after the above fix, set it to zero
    if math.isnan(min_mnll):
        min_mnll = 0.0

    if math.isinf(min_mnll):
        min_mnll = 0.0

    # mnll of profile with uniform profile
    max_mnll = mnll(profile, probs=uniform_profile)

    return (min_mnll, max_mnll)

#https://github.com/kundajelab/basepairmodels/blob/cf8e346e9df1bad9e55bd459041976b41207e6e5/basepairmodels/cli/fastpredict.py#L131
def jsd_min_max_bounds(profile):
    """
        Min Max bounds for the jsd metric
        
        Args:
            profile (numpy.ndarray): true profile 
            
        Returns:
            tuple: (min, max) bounds values
    """
    
    # uniform distribution profile
    uniform_profile = np.ones(len(profile)) * (1.0 / len(profile))

    # profile as probabilities
    profile_prob = profile / np.sum(profile)

    # jsd of profile with uniform profile
    max_jsd = jensenshannon(profile_prob, uniform_profile)

    # jsd of profile with itself (upper bound)
    min_jsd = 0.0

    return (min_jsd, max_jsd)
