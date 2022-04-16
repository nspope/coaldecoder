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

Multithreaded version
TODO: this still uses GIL, and I don't think reticulate works with multiprocessing
"""

import numpy as np
import scipy as spy
import tskit
from math import factorial as fac
from itertools import combinations_with_replacement as combn_replace
from itertools import combinations as combn
from collections import namedtuple
import concurrent.futures

#-------------------- streamlined rate calculation without edge diffs --------------------#

NodeWeightTable = namedtuple('NodeWeightTable', ['id', 'time', 'weights'])

class ObservedTrioRatesMultithreaded:
    """
    Class that calculates rates of first and second trio coalescences in a tree
    sequence and optionally block-bootstraps these counts
    """

    def __init__ (self, ts, sample_sets, time_breaks, bootstrap_replicates=0, bootstrap_blocks=1, random_seed=None, mask=None, threads=1):

        self.prefix = "[CalculateTrioRates] "

        # check inputs
        if isinstance(ts, str):
            ts = tskit.load(ts)
        else:
            assert isinstance(ts, tskit.TreeSequence)

        assert isinstance(sample_sets, dict)
        self.population_names = np.sort(list(sample_sets))
        self.sample_sets = []
        for name in self.population_names:
            s = sample_sets[name]
            assert all([isinstance(i, int) for i in s])
            assert all([i in ts.samples() for i in s])
            self.sample_sets.append(s)
        self.num_populations = len(self.sample_sets)

        #TODO:
        # what should time breaks instance be?

        bootstrap_blocks = int(bootstrap_blocks)
        assert bootstrap_blocks > 0

        bootstrap_replicates = int(bootstrap_replicates)
        assert bootstrap_replicates >= 0

        if random_seed is not None:
            random_seed = int(random_seed)
            assert random_seed >= 0

        if mask is not None:
            assert isinstance(mask, np.ndarray)
            print(self.prefix + "Removing masked intervals from tree sequence")
            ts.delete_intervals(mask)

        threads = int(threads)
        assert threads > 0

        self.sequence_length = ts.sequence_length
        self.num_trees = ts.num_trees

        # make pair indices
        self._pair_index = np.zeros((self.num_populations, self.num_populations), "int32")
        for i, u in enumerate(combn_replace(np.arange(self.num_populations), 2)):
            self._pair_index[u[0], u[1]] = i
        #TODO: lambda can't be pickled, can I use regular function?
        self._pair_linear_index = lambda i, j : self._pair_index[i, j]

        # check that sample times are all the same within a sample set
        sample_times = [list({ts.get_time(i) for i in s}) for s in self.sample_sets]
        assert all([len(t) == 1 for t in sample_times])
        self.population_times = [t[0] for t in sample_times]

        # check that time breaks include sample times
        assert all([i in time_breaks for i in self.population_times])
        self.time_breaks = time_breaks

        # TODO: allow user to provide start/end of genomic interval
        # calculate tree and block span
        tree_idx = np.arange(ts.first().index, ts.last().index + 1)
        bootstrap_blocks = min(bootstrap_blocks, len(tree_idx))
        tree_block = np.floor_divide(bootstrap_blocks * tree_idx, len(tree_idx))
        tree_span = np.array([
            tree.interval.right - tree.interval.left for tree in ts.trees()
        ])
        tree_span /= np.sum(tree_span)
        block_span = np.array([
            sum(tree_span[np.where(tree_block == i)]) for i in range(bootstrap_blocks)
        ])

        # calculate weights per block
        # TODO: For lots of pops the size of the weights table will be huge. It could be a sparse matrix instead
        # TODO: use edge diffs
        num_weights = 2*int(self.num_populations*self.num_populations*(self.num_populations+1)/2)
        num_time_breaks = len(self.time_breaks)
        block_y = np.zeros((num_weights, num_time_breaks, bootstrap_blocks))
        block_n = np.zeros((num_weights, num_time_breaks, bootstrap_blocks))
        
        def _weights_per_block (_block):
            _min_tree = tree_idx[np.searchsorted(tree_block, _block, side='left')]
            _max_tree = tree_idx[np.searchsorted(tree_block, _block, side='right') - 1]
            _tree = ts.at_index(_min_tree)
            _y = np.zeros((num_weights, num_time_breaks))
            _n = np.zeros((num_weights, num_time_breaks))
            while _tree.index <= _max_tree:
                _weights = self._trio_counts_for_tree(_tree)
                _time_bin = np.searchsorted(self.time_breaks, _weights.time, side='right') - 1
                _norm_weights = _weights.weights * tree_span[i]
                for j in range(len(_weights.id)):
                    #block_y[:, _time_bin[j], _block] += _norm_weights[:, j]
                    _y[:, _time_bin[j]] += _norm_weights[:, j]
                    for k in range(_time_bin[j] + 1):
                        #block_n[:, k, _block] += _norm_weights[:, j]
                        _n[:, k] += _norm_weights[:, j]
                if _tree.index == _max_tree: break
                _tree.next()

        #for i in range(bootstrap_blocks): _weights_per_block(i)
        with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
            executor.map(_weights_per_block, range(bootstrap_blocks))

        block_n = self._share_denominator_across_initial_states(block_n)

        #TODO: move this to rates()
        self.y = np.zeros((num_weights, num_time_breaks))
        self.n = np.zeros((num_weights, num_time_breaks))
        for block in range(bootstrap_blocks):
            self.y[:,:] += block_y[:,:,block]
            self.n[:,:] += block_n[:,:,block]

        #TODO: move this to bootstrap_rates()
        # bootstrap
        rng = np.random.default_rng(random_seed)
        self.y_boot = np.zeros((num_weights, num_time_breaks, bootstrap_replicates))
        self.n_boot = np.zeros((num_weights, num_time_breaks, bootstrap_replicates))
        for rep in range(bootstrap_replicates):
            block_multiplier = rng.multinomial(
                bootstrap_blocks, [1.0 / bootstrap_blocks] * bootstrap_blocks
            )
            block_weights = block_span * block_multiplier
            block_weights /= np.sum(block_weights)
            for block in range(bootstrap_blocks):
                self.y_boot[:,:,rep] += block_weights[block]/block_span[block] * block_y[:,:,block]
                self.n_boot[:,:,rep] += block_weights[block]/block_span[block] * block_n[:,:,block]

    def _trio_counts_for_tree (self, tree):
        """
        Count numbers of first and second coalescent events in a tree and return in a table.
        The "weights" field of this table has nodes for columns and emissions for rows.
        """

        sample_sets = self.sample_sets

        # ignore isolated samples
        nodes = np.sort([i for i in tree.nodes() if not tree.is_isolated(i)])
        samples = np.sort([i for i in tree.samples() if not tree.is_isolated(i)])

        P = len(sample_sets)
        N = len(nodes)
        S = [np.intersect1d(samples, s).size for s in sample_sets]

        # number of descendant/non-descendant samples in each population
        n_desc = np.zeros((N, P))
        n_outg = np.zeros((N, P))
        nz_desc = []
        nz_outg = []
        time = np.zeros((N))
        for i, node in enumerate(nodes):
            si = [j for j in tree.samples(node) if not tree.is_isolated(j)]
            n_desc[i, :] = [np.intersect1d(si, s).size for s in sample_sets]
            n_outg[i, :] = [s - n for s, n in zip(S, n_desc[i, :])]
            nz_desc += [np.where(n_desc[i, :] > 0)[0]]
            nz_outg += [np.where(n_outg[i, :] > 0)[0]]
            time[i] = tree.time(node)

        # number of sample pairs that coalesce at node
        pairs = int(P*(P+1)/2)
        n_pairs_at = np.zeros((N, pairs))
        nz_pairs_at = []
        for i, node in enumerate(nodes):
            ch = np.searchsorted(nodes, tree.children(node))
            for a, b in combn_replace(nz_desc[i], 2):
                j = self._pair_linear_index(a, b)
                for u, v in combn(ch, 2):
                    n_pairs_at[i, j] += (n_desc[u,a]*n_desc[v,b] + n_desc[v,a]*n_desc[u,b]) / (1+int(a==b))
            nz_pairs_at += [np.where(n_pairs_at[i, :] > 0)[0]]

        # number of sample pairs that coalesce at or above node
        n_pairs_above = np.zeros((N, pairs))
        nz_pairs_above = []
        for i, node in enumerate(nodes):
            ch = np.searchsorted(nodes, tree.children(node))
            for a, b in combn_replace(nz_desc[i], 2):
                j = self._pair_linear_index(a, b)
                if a == b:
                    n_pairs_above[i, j] += n_desc[i,a]*(n_desc[i,a]-1)/2
                else:
                    n_pairs_above[i, j] += n_desc[i,a]*n_desc[i,b]
            nz_pairs_above += [np.where(n_pairs_above[i, :] > 0)[0]]

        # number of first/second trio coalescence events
        trios = int(P*P*(P+1)/2)
        n_first_trios = np.zeros((N, trios))
        n_second_trios = np.zeros((N, trios))
        for i, node in enumerate(nodes):
            ch = np.searchsorted(nodes, tree.children(node))

            # first coalescence events
            for j in nz_pairs_at[i]:
                for c in nz_outg[i]:
                    k = j * P + c
                    n_first_trios[i, k] += n_pairs_at[i, j] * n_outg[i, c]

            # second coalescence events
            for u, v in combn(ch, 2): 
                for j in nz_pairs_above[u]:
                    for c in nz_desc[v]:
                        k = j * P + c
                        n_second_trios[i, k] += n_pairs_above[u, j] * n_desc[v, c]
                for j in nz_pairs_above[v]:
                    for c in nz_desc[u]:
                        k = j * P + c
                        n_second_trios[i, k] += n_pairs_above[v, j] * n_desc[u, c]

            for u, v, w in combn(ch, 3): 
                #if node is ternary or greater:
                #do we need to consider star-coalescences?
                #in this case can't say which pair coalesced first so ignore
                pass

        weights = np.concatenate([n_first_trios.T, n_second_trios.T], 0)
        
        return NodeWeightTable(nodes, time, weights)

    def _share_denominator_across_initial_states (self, n):
        """
        Sum denominator of rate for emissions with the same initial state, then
        map back to original array
        """

        denominator_labels = self.denominator_labels()
        unique_denominator_labels = np.sort(np.unique(denominator_labels))
        denominator_labels_idx = np.searchsorted(unique_denominator_labels, denominator_labels, side='right') - 1
        n_pool = np.zeros((len(unique_denominator_labels), n.shape[1], n.shape[2]))
        for i, j in enumerate(denominator_labels_idx):
            n_pool[j, :, :] += n[i, :, :]
        n_share = np.zeros(n.shape)
        for i, j in enumerate(denominator_labels_idx):
            n_share[i, :, :] = n_pool[j, :, :]
        return n_share

    def emission_states (self):
        trio_labels = []
        for t in ['t1::', 't2::']:
            for a, b in combn_replace(range(self.num_populations), 2):
                pair_label = t + '((' + self.population_names[a] + ',' + self.population_names[b] + '),'
                for c in range(self.num_populations):
                    trio_labels += [pair_label + self.population_names[c] + ")"]
        return trio_labels

    def denominator_labels (self):
        denominator_labels = []
        for t in ['t1::', 't2::']:
            for a, b in combn_replace(range(self.num_populations), 2):
                pair = [a, b]
                for c in range(self.num_populations):
                    u = np.sort(pair + [c]) #if t == 't1::' else pair + [c]
                    denominator_labels.append(t + "{" + 
                            self.population_names[u[0]] + "," + 
                            self.population_names[u[1]] + "," + 
                            self.population_names[u[2]] + "}")
        return denominator_labels

    def epochs (self):
        return ["[" + str(self.time_breaks[i]) + "," + str(self.time_breaks[i+1]) + ")" 
                for i in range(len(self.time_breaks)-1)]

    def join (self, rhs):
        assert isinstance (rhs, ObservedTrioRates)
        assert all([i == j for i, j in zip(self.epochs(), rhs.epochs())])
        assert all([i == j for i, j in zip(self.emission_states(), rhs.emission_states())])

        total_length = self.sequence_length + rhs.sequence_length
        lhs_weight = self.sequence_length / total_length
        rhs_weight = rhs.sequence_length / total_length

        self.y = self.y*lhs_weight + rhs.y*rhs_weight
        self.n = self.n*lhs_weight + rhs.n*rhs_weight
        self.y_boot = self.y_boot*lhs_weight + rhs.y_boot*rhs_weight
        self.n_boot = self.n_boot*lhs_weight + rhs.n_boot*rhs_weight
        self.sequence_length = total_length
        self.num_trees += rhs.num_trees

    def rates (self):
        rates = self.y[:,:-1] / self.n[:,:-1]
        for i in range(rates.shape[1]):
            rates[:,i] /= (self.time_breaks[i+1] - self.time_breaks[i])
        return rates

    def bootstrapped_rates (self):
        rates_boot = self.y_boot[:,:-1,:] / self.n_boot[:,:-1,:]
        for i in range(rates_boot.shape[1]):
            rates_boot[:,i,:] /= (self.time_breaks[i+1] - self.time_breaks[i])
        return rates_boot
