#
# MIT License
#
# Copyright (c) 2021-2022 Nathaniel S. Pope
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Helper functions to extract trio coalescence rates from tree sequences.
"""

import numpy as np
import tskit
import copy
import itertools

def trio_first_coalescence_rates(
    ts_paths,
    sample_sets,
    time_breaks,
    bootstrap_replicates = 0,
    bootstrap_blocks = 1,
    random_seed = None,
    ):

    time_breaks = np.append(np.array(time_breaks, dtype="float"), 0)
    time_breaks = np.sort(np.unique(time_breaks))

    bootstrap_replicates = int(bootstrap_replicates)
    bootstrap_blocks = int(bootstrap_blocks)
    if random_seed is not None:
        random_seed = int(random_seed)

    labels = []
    for a,b in itertools.combinations_with_replacement(range(len(sample_sets)), 2):
        for c in range(len(sample_sets)):
            labels += [[a,b,c]]

    out = []
    for path in ts_paths:
        # read tree sequence
        ts = tskit.load(path)
        distr = CoalescenceTimeDistribution(
            ts,
            weight_func='trio_first_coalescence_events', 
            sample_sets=sample_sets,
            blocks_per_window=bootstrap_blocks,
            span_normalise=True
        )
        n = np.nan_to_num(distr.num_uncoalesced(np.array([0.0])))
        n = n[:,0,0]
        y = np.nan_to_num(distr.num_coalesced(time_breaks))
        y = y[:,1:,0] - y[:,:-1,0]
        y_boot = np.zeros((y.shape[0], y.shape[1], bootstrap_replicates))
        n_boot = np.zeros((n.shape[0], bootstrap_replicates))
        if bootstrap_replicates > 0:
            generator = distr.block_bootstrap(
                    random_seed=random_seed, 
                    num_replicates=bootstrap_replicates)
            for distr_boot in generator:
                n_tmp = np.nan_to_num(distr_boot.num_uncoalesced(np.array([0.0])))
                n_boot[:,i] = n_tmp[:,0,0]
                y_tmp = np.nan_to_num(distr_boot.num_coalesced(time_breaks))
                y_boot[:,:,i] = y_tmp[:,1:,0] - y_tmp[:,:-1,0]
        out += [{
            "file" : path, "y" : y, "n" : n, 
            "y_boot" : y_boot, "n_boot" : n_boot, 
            "size" : ts.sequence_length, "trees" : ts.num_trees, 
            "labels" : np.array(labels),
        }]
    return out
