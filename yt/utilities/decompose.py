"""
Automagical cartesian domain decomposition.

Author: Kacper Kowalik <xarthisius.kk@gmail.com>
Affiliation: CA UMK
Author: Artur Gawryszczak <gawrysz@gmail.com>
Affiliation: PCSS
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Kacper Kowalik. All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

SIEVE_PRIMES = \
    lambda l: l and l[:1] + SIEVE_PRIMES([n for n in l if n % l[0]])


def decompose_to_primes(max_prime):
    """ Decompose number into the primes """
    for prime in SIEVE_PRIMES(range(2, max_prime)):
        if prime * prime > max_prime:
            break
        while max_prime % prime == 0:
            yield prime
            max_prime /= prime
    if max_prime > 1:
        yield max_prime


def decompose_array(arr, psize, bbox):
    """ Calculate list of product(psize) subarrays of arr, along with their
        left and right edges
    """
    grid_left_edges = np.empty([np.product(psize), 3], dtype=np.float64)
    grid_right_edges = np.empty([np.product(psize), 3], dtype=np.float64)
    n_d = arr.shape
    d_s = (bbox[:, 1] - bbox[:, 0]) / n_d
    dist = np.mgrid[bbox[0, 0]:bbox[0, 1]:d_s[0],
                    bbox[1, 0]:bbox[1, 1]:d_s[1],
                    bbox[2, 0]:bbox[2, 1]:d_s[2]]
    for i in range(3):
        xyz = split_array(dist[i], psize)
        for j in range(np.product(psize)):
            grid_left_edges[j, i] = xyz[j][0, 0, 0]
            grid_right_edges[j, i] = xyz[j][-1, -1, -1] + d_s[i]
        del xyz
    del dist
    patches = split_array(arr, psize)
    return grid_left_edges, grid_right_edges, patches


def evaluate_domain_decomposition(n_d, pieces, ldom):
    """ Evaluate longest to shortest edge ratio
        BEWARE: lot's of magic here """
    eff_dim = (n_d > 1).sum()
    ideal_bsize = eff_dim * (pieces * np.product(n_d) ** (eff_dim - 1)
                             ) ** (1.0 / eff_dim)
    mask = np.where(n_d > 1)
    nd = np.array(n_d, dtype=np.float64)[mask]
    bsize = int(np.sum(ldom[mask] / nd * np.product(nd)))
    load_balance = float(np.product(n_d)) / \
        (float(pieces) * np.product((n_d - 1) / ldom + 1))

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
        [(prime, temp[prime], (temp[prime] + 1) * (temp[prime] + 2) / 2)
         for prime in np.unique(factors)]
    )


def get_psize(n_d, pieces):
    """ Calculate the best division of array into px*py*pz subarrays.
        The goal is to minimize the ratio of longest to shortest edge
        to minimize the amount of inter-process communication.
    """
    fac = factorize_number(pieces)
    nfactors = len(fac[:, 2])
    best = 0.0
    while np.all(fac[:, 2] > 0):
        ldom = np.ones(3, dtype=np.int)
        for nfac in range(nfactors):
            i = int(np.sqrt(0.25 + 2 * (fac[nfac, 2] - 1)) - 0.5)
            k = fac[nfac, 2] - int(1 + i * (i + 1) / 2)
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


def split_array(tab, psize):
    """ Split array into px*py*pz subarrays. """
    n_d = np.array(tab.shape, dtype=np.int64)
    slices = []
    for i in range(psize[0]):
        for j in range(psize[1]):
            for k in range(psize[2]):
                pc = np.array((i, j, k), dtype=np.int64)
                le = n_d * pc / psize
                re = n_d * (pc + np.ones(3, dtype=np.int64)) / psize
                slices.append(np.s_[le[0]:re[0], le[1]:re[1], le[2]:re[2]])
    return [tab[sl] for sl in slices]
