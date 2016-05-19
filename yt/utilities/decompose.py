"""
Automagical cartesian domain decomposition.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division

import numpy as np

SIEVE_PRIMES = \
    lambda l: l and l[:1] + SIEVE_PRIMES([n for n in l if n % l[0]])


def decompose_to_primes(max_prime):
    """ Decompose number into the primes """
    for prime in SIEVE_PRIMES(list(range(2, max_prime))):
        if prime * prime > max_prime:
            break
        while max_prime % prime == 0:
            yield prime
            max_prime //= prime
    if max_prime > 1:
        yield max_prime

def decompose_array(shape, psize, bbox):
    """ Calculate list of product(psize) subarrays of arr, along with their
        left and right edges
    """
    return split_array(bbox[:, 0], bbox[:, 1], shape, psize)


def evaluate_domain_decomposition(n_d, pieces, ldom):
    """ Evaluate longest to shortest edge ratio
        BEWARE: lot's of magic here """
    eff_dim = (n_d > 1).sum()
    exp = float(eff_dim - 1) / float(eff_dim)
    ideal_bsize = eff_dim * pieces ** (1.0 / eff_dim) * np.product(n_d) ** exp
    mask = np.where(n_d > 1)
    nd_arr = np.array(n_d, dtype=np.float64)[mask]
    bsize = int(np.sum(ldom[mask] / nd_arr * np.product(nd_arr)))
    load_balance = float(np.product(n_d)) / \
        (float(pieces) * np.product((n_d - 1) // ldom + 1))

    # 0.25 is magic number
    quality = load_balance / (1 + 0.25 * (bsize / ideal_bsize - 1.0))
    # \todo add a factor that estimates lower cost when x-direction is
    # not chopped too much
    # \deprecated estimate these magic numbers
    quality *= (1. - (0.001 * ldom[0] + 0.0001 * ldom[1]) / pieces)
    if np.any(ldom > n_d):
        quality = 0

    return quality


def factorize_number(pieces):
    """ Return array consiting of prime, its power and number of different
        decompositions in three dimensions for this prime
    """
    factors = [factor for factor in decompose_to_primes(pieces)]
    temp = np.bincount(factors)
    return np.array(
        [(prime, temp[prime], (temp[prime] + 1) * (temp[prime] + 2) // 2)
         for prime in np.unique(factors)]).astype(np.int64)


def get_psize(n_d, pieces):
    """ Calculate the best division of array into px*py*pz subarrays.
        The goal is to minimize the ratio of longest to shortest edge
        to minimize the amount of inter-process communication.
    """
    fac = factorize_number(pieces)
    nfactors = len(fac[:, 2])
    best = 0.0
    p_size = np.ones(3, dtype=np.int64)
    if pieces == 1:
        return p_size

    while np.all(fac[:, 2] > 0):
        ldom = np.ones(3, dtype=np.int64)
        for nfac in range(nfactors):
            i = int(np.sqrt(0.25 + 2 * (fac[nfac, 2] - 1)) - 0.5)
            k = fac[nfac, 2] - (1 + i * (i + 1) // 2)
            i = fac[nfac, 1] - i
            j = fac[nfac, 1] - (i + k)
            ldom *= fac[nfac, 0] ** np.array([i, j, k])

        quality = evaluate_domain_decomposition(n_d, pieces, ldom)
        if quality > best:
            best = quality
            p_size = ldom
        # search for next unique combination
        for j in range(nfactors):
            if fac[j, 2] > 1:
                fac[j, 2] -= 1
                break
            else:
                if (j < nfactors - 1):
                    fac[j, 2] = int((fac[j, 1] + 1) * (fac[j, 1] + 2) / 2)
                else:
                    fac[:, 2] = 0  # no more combinations to try

    return p_size

def split_array(gle, gre, shape, psize):
    """ Split array into px*py*pz subarrays. """
    n_d = np.array(shape, dtype=np.int64)
    dds = (gre-gle)/shape
    left_edges = []
    right_edges = []
    shapes = []
    slices = []
    for i in range(psize[0]):
        for j in range(psize[1]):
            for k in range(psize[2]):
                piece = np.array((i, j, k), dtype=np.int64)
                lei = n_d * piece // psize
                rei = n_d * (piece + np.ones(3, dtype=np.int64)) // psize
                lle = gle + lei*dds
                lre = gle + rei*dds
                left_edges.append(lle)
                right_edges.append(lre)
                shapes.append(rei-lei)
                slices.append(np.s_[lei[0]:rei[0], lei[1]:
                                    rei[1], lei[2]:rei[2]])

    return left_edges, right_edges, shapes, slices
