"""Dinucleotide-preserving shuffles used as DeepSHAP references.

The algorithm is a verbatim port of deeplift 0.6.13 ``deeplift/dinuc_shuffle.py``
(https://github.com/kundajelab/deeplift), so a given ``numpy.random.RandomState`` produces exactly the
references deeplift produces. Changes: ``rng`` is explicit (``None`` still means a fresh unseeded RandomState,
as in deeplift), ``tostring`` -> ``tobytes`` for NumPy 2, and an optional numba version of the Eulerian walk
(the permutations still come from the same RandomState calls, in the same order, so results are identical).

deeplift is distributed under the MIT License:

    MIT License

    Copyright (c) 2018 Kundaje Lab

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
"""
import zlib

import numpy as np

try:
    import numba
except ImportError:  # pragma: no cover - numba is a regular dependency (modisco needs it)
    numba = None

NUM_SHUFFLES = 20


def string_to_char_array(seq):
    """
    Converts an ASCII string to a NumPy array of byte-long ASCII codes.
    e.g. "ACGT" becomes [65, 67, 71, 84].
    """
    return np.frombuffer(bytearray(seq, "utf8"), dtype=np.int8)


def char_array_to_string(arr):
    """
    Converts a NumPy array of byte-long ASCII codes into an ASCII string.
    e.g. [65, 67, 71, 84] becomes "ACGT".
    """
    return arr.tobytes().decode("ascii")


def one_hot_to_tokens(one_hot):
    """
    Converts an L x D one-hot encoding into an L-vector of integers in the range
    [0, D], where the token D is used when the one-hot encoding is all 0. This
    assumes that the one-hot encoding is well-formed, with at most one 1 in each
    column (and 0s elsewhere).
    """
    tokens = np.tile(one_hot.shape[1], one_hot.shape[0])  # Vector of all D
    seq_inds, dim_inds = np.where(one_hot)
    tokens[seq_inds] = dim_inds
    return tokens


def tokens_to_one_hot(tokens, one_hot_dim):
    """
    Converts an L-vector of integers in the range [0, D] to an L x D one-hot
    encoding. The value `D` must be provided as `one_hot_dim`. A token of D
    means the one-hot encoding is all 0s.
    """
    identity = np.identity(one_hot_dim + 1)[:, :-1]  # Last row is all 0s
    return identity[tokens]


def _walk_python(tokens, shuf_next_inds):
    counters = [0] * len(shuf_next_inds)

    # Build the resulting array
    ind = 0
    result = np.empty_like(tokens)
    result[0] = tokens[ind]
    for j in range(1, len(tokens)):
        t = tokens[ind]
        ind = shuf_next_inds[t][counters[t]]
        counters[t] += 1
        result[j] = tokens[ind]
    return result


def _walk_kernel_py(tokens, flat_next_inds, offsets, result):
    counters = np.zeros(offsets.shape[0] - 1, dtype=np.int64)
    ind = 0
    result[0] = tokens[ind]
    for j in range(1, tokens.shape[0]):
        t = tokens[ind]
        k = offsets[t] + counters[t]
        if k >= offsets[t + 1]:
            return False
        ind = flat_next_inds[k]
        counters[t] += 1
        result[j] = tokens[ind]
    return True


# numba's on-disk cache needs a writable cache directory; without one, njit(cache=True) raises at decoration time
_walk_kernel = None
if numba is not None:
    for _cache in (True, False):
        try:
            _walk_kernel = numba.njit(cache=_cache, nogil=True)(_walk_kernel_py)
            break
        except (RuntimeError, OSError):
            continue


def _walk_numba(tokens, shuf_next_inds):
    tokens = np.ascontiguousarray(tokens, dtype=np.int64)
    sizes = np.array([len(inds) for inds in shuf_next_inds], dtype=np.int64)
    offsets = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=offsets[1:])
    flat = np.concatenate(shuf_next_inds).astype(np.int64)
    result = np.empty_like(tokens)
    if not _walk_kernel(tokens, flat, offsets, result):
        raise IndexError("dinucleotide shuffle walk ran out of successors")
    return result


def dinuc_shuffle(seq, num_shufs=None, rng=None, use_numba=True):
    """
    Creates shuffles of the given sequence, in which dinucleotide frequencies
    are preserved.
    Arguments:
        `seq`: either a string of length L, or an L x D NumPy array of one-hot
            encodings
        `num_shufs`: the number of shuffles to create, N; if unspecified, only
            one shuffle will be created
        `rng`: a NumPy RandomState object, to use for performing shuffles
        `use_numba`: run the walk with numba when available (identical results)
    If `seq` is a string, returns a list of N strings of length L, each one
    being a shuffled version of `seq`. If `seq` is a 2D NumPy array, then the
    result is an N x L x D NumPy array of shuffled versions of `seq`, also
    one-hot encoded. If `num_shufs` is not specified, then the first dimension
    of N will not be present (i.e. a single string will be returned, or an L x D
    array).
    """
    if type(seq) is str:
        arr = string_to_char_array(seq)
    elif type(seq) is np.ndarray and len(seq.shape) == 2:
        seq_len, one_hot_dim = seq.shape
        arr = one_hot_to_tokens(seq)
    else:
        raise ValueError("Expected string or one-hot encoded array")

    if rng is None:
        rng = np.random.RandomState()

    walk = _walk_numba if (use_numba and _walk_kernel is not None) else _walk_python

    # Get the set of all characters, and a mapping of which positions have which
    # characters; use `tokens`, which are integer representations of the
    # original characters
    chars, tokens = np.unique(arr, return_inverse=True)
    tokens = tokens.reshape(-1)

    # For each token, get a list of indices of all the tokens that come after it
    shuf_next_inds = []
    for t in range(len(chars)):
        mask = tokens[:-1] == t  # Excluding last char
        inds = np.where(mask)[0]
        shuf_next_inds.append(inds + 1)  # Add 1 for next token

    if type(seq) is str:
        all_results = []
    else:
        all_results = np.empty(
            (num_shufs if num_shufs else 1, seq_len, one_hot_dim),
            dtype=seq.dtype
        )

    for i in range(num_shufs if num_shufs else 1):
        # Shuffle the next indices
        for t in range(len(chars)):
            inds = np.arange(len(shuf_next_inds[t]))
            inds[:-1] = rng.permutation(len(inds) - 1)  # Keep last index same
            shuf_next_inds[t] = shuf_next_inds[t][inds]

        result = walk(tokens, shuf_next_inds)

        if type(seq) is str:
            all_results.append(char_array_to_string(chars[result]))
        else:
            all_results[i] = tokens_to_one_hot(chars[result], one_hot_dim)
    return all_results if num_shufs else all_results[0]


def reference_seed(onehot, seed=None):
    """RandomState seed for one sequence: crc32 of its int8 one-hot bytes, xor `seed` when given.

    Seeding from the content makes a sequence's references independent of its position in the input, of the
    subsample and of the batching. `seed=None` is the convention of the legacy parity goldens
    (tests/goldens/make_goldens.py content_seed).
    """
    crc = zlib.crc32(np.ascontiguousarray(onehot, dtype=np.int8).tobytes()) & 0xFFFFFFFF
    if seed is None:
        return crc
    return (crc ^ int(seed)) & 0xFFFFFFFF


def make_references(seqs, num_shuffles=NUM_SHUFFLES, seed=1234, use_numba=True):
    """`num_shuffles` dinucleotide shuffles per sequence: (N, L, 4) one-hot -> (N, num_shuffles, L, 4) int8.

    Each sequence gets one dinuc_shuffle call with num_shufs=num_shuffles and its own
    RandomState(reference_seed(sequence, seed)).
    """
    seqs = np.asarray(seqs)
    if seqs.ndim != 3:
        raise ValueError("expected one-hot sequences of shape (N, L, 4), got {}".format(seqs.shape))
    refs = np.empty((seqs.shape[0], num_shuffles) + seqs.shape[1:], dtype=np.int8)
    for i in range(seqs.shape[0]):
        s = np.ascontiguousarray(seqs[i], dtype=np.int8)
        rng = np.random.RandomState(reference_seed(s, seed))
        refs[i] = dinuc_shuffle(s, num_shufs=num_shuffles, rng=rng, use_numba=use_numba)
    return refs
