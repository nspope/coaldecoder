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
import scipy as spy
import tskit
from math import factorial as fac
from itertools import combinations_with_replacement as combn_replace
from itertools import combinations as combn
from collections import namedtuple

#-------------------- streamlined rate calculation without edge diffs --------------------#

NodeWeightTable = namedtuple('NodeWeightTable', ['id', 'time', 'min_time', 'weights'])

class ObservedTrioRates:
    """
    Class that calculates rates of first and second trio coalescences in a tree
    sequence and optionally block-bootstraps these counts
    """

    def __init__ (self, ts, sample_sets, time_breaks, bootstrap_blocks=1, trees_per_block=None, mask=None, threads=1):

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

        if trees_per_block is not None:
            assert trees_per_block > 0
            bootstrap_blocks = ts.num_trees / trees_per_block

        bootstrap_blocks = int(bootstrap_blocks)
        assert bootstrap_blocks > 0

        #TODO:
        # alternatively specify trees_per_block

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
        #self._pair_linear_index = lambda i, j : self._pair_index[i, j]

        # check that sample times are all the same within a sample set
        sample_times = [list({ts.get_time(i) for i in s}) for s in self.sample_sets]
        assert all([len(t) == 1 for t in sample_times])
        self.population_times = [t[0] for t in sample_times]

        # check that time breaks include sample times
        assert all([i in time_breaks for i in self.population_times])
        self.time_breaks = time_breaks

        # TODO: allow user to provide start/end of genomic interval
        # calculate tree and block span
        tree_exists = np.array([tree.num_edges > 0 for tree in ts.trees()])
        tree_idx = np.array([-1 for _ in range(ts.num_trees)])
        tree_idx[tree_exists] = np.arange(np.sum(tree_exists))
        bootstrap_blocks = min(bootstrap_blocks, np.sum(tree_exists))
        tree_block = np.floor_divide(bootstrap_blocks * tree_idx, np.sum(tree_exists))
        tree_span = np.array([
            tree.interval.right - tree.interval.left for tree in ts.trees()
        ])
        tree_span *= tree_exists
        tree_span /= np.sum(tree_span)
        block_span = np.array([
            sum(tree_span[np.where(tree_block == i)]) for i in range(bootstrap_blocks)
        ])

        # TODO: For lots of pops the size of the weights table will be huge. It could be a sparse matrix instead
        # TODO: parallelize over blocks
        # TODO: use edge diffs
        # calculate weights per block
        num_weights = 2*int(self.num_populations*self.num_populations*(self.num_populations+1)/2)
        num_time_breaks = len(self.time_breaks)
        block_y = np.zeros((num_weights, num_time_breaks, bootstrap_blocks))
        block_n = np.zeros((num_weights, num_time_breaks, bootstrap_blocks))
        for i, tree in enumerate(ts.trees()):
            if tree_exists[i]:
                block = tree_block[i]
                weights = self._trio_counts_for_tree(tree)
                time_bin = np.searchsorted(self.time_breaks, weights.time, side='right') - 1
                norm_weights = weights.weights * tree_span[i]
                for j in range(len(weights.id)):
                    block_y[:, time_bin[j], block] += norm_weights[:, j]
                    for k in range(time_bin[j] + 1):
                        block_n[:, k, block] += norm_weights[:, j]

        # TODO: initialize all members here, not piecemeal
        self.n = self._share_denominator_across_initial_states(block_n)
        self.y = block_y
        self.block_span = block_span
        self.bootstrap_blocks = bootstrap_blocks
        self.num_weights = num_weights
        self.num_time_breaks = num_time_breaks

    def _pair_linear_index (self, i, j):
        return self._pair_index[i, j]

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
        min_time = np.zeros((N))
        time = np.zeros((N))
        for i, node in enumerate(nodes):
            si = [j for j in tree.samples(node) if not tree.is_isolated(j)]
            n_desc[i, :] = [np.intersect1d(si, s).size for s in sample_sets]
            n_outg[i, :] = [s - n for s, n in zip(S, n_desc[i, :])]
            nz_desc += [np.where(n_desc[i, :] > 0)[0]]
            nz_outg += [np.where(n_outg[i, :] > 0)[0]]
            min_time[i] = max([tree.time(j) for j in si])
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
        
        return NodeWeightTable(nodes, time, min_time, weights)

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
                    u = np.sort(pair + [c])
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

        self.y = np.concatenate((self.y*lhs_weight, rhs.y*rhs_weight), axis=2)
        self.n = np.concatenate((self.n*lhs_weight, rhs.n*rhs_weight), axis=2)
        self.block_span = np.concatenate((self.block_span*lhs_weight, rhs.block_span*rhs_weight))

        self.bootstrap_blocks += rhs.bootstrap_blocks
        self.sequence_length = total_length
        self.num_trees += rhs.num_trees

    def numerator (self, block_weights = None):
        if block_weights is None:
            block_weights = np.ones((self.bootstrap_blocks))
            block_weights /= np.sum(block_weights)
        else:
            assert len(block_weights) == self.bootstrap_blocks
            assert np.all(block_weights >= 0.0)
        y = np.zeros((self.num_weights, self.num_time_breaks))
        for block in range(self.bootstrap_blocks):
            y[:,:] += block_weights[block] * self.y[:,:,block]
        return y[:,:-1]

    def denominator (self, block_weights = None):
        if block_weights is None:
            block_weights = np.ones((self.bootstrap_blocks))
            block_weights /= np.sum(block_weights)
        else:
            assert len(block_weights) == self.bootstrap_blocks
            assert np.all(block_weights >= 0.0)
        n = np.zeros((self.num_weights, self.num_time_breaks))
        for block in range(self.bootstrap_blocks):
            n[:,:] += block_weights[block] * self.n[:,:,block]
        return n[:,:-1]

    def rates (self, block_weights = None):
        rates = self.numerator(block_weights) / self.denominator(block_weights)
        for i in range(self.num_time_breaks - 1):
            rates[:,i] /= (self.time_breaks[i+1] - self.time_breaks[i])
        return rates

    def block_bootstrap (self, num_replicates, random_seed = None):
        num_replicates = int(num_replicates)
        assert num_replicates > 0

        if random_seed is not None:
            random_seed = int(random_seed)
            assert random_seed >= 0

        rng = np.random.default_rng(random_seed)
        for rep in range(num_replicates):
            block_multiplier = rng.multinomial(
                self.bootstrap_blocks, [1.0 / self.bootstrap_blocks] * self.bootstrap_blocks
            )
            block_weights = self.block_span * block_multiplier
            block_weights /= np.sum(block_weights)
            rates = self.rates(block_weights / self.block_span)
            yield rates

    def std_dev (self, num_replicates, random_seed = None):
        num_replicates = int(num_replicates)
        assert num_replicates > 1

        n = 0.0
        mean = np.zeros((self.num_weights, self.num_time_breaks-1))
        sumsq = np.zeros((self.num_weights, self.num_time_breaks-1))
        for x in self.block_bootstrap(num_replicates, random_seed):
            n += 1.0
            delta = x - mean
            mean += delta / n
            delta *= x - mean
            sumsq += delta
        return np.sqrt(sumsq / (n - 1.0))

    def bootstrapped_rates (self, num_replicates, random_seed = None):
        num_replicates = int(num_replicates)
        assert num_replicates > 0

        rates = np.zeros((self.num_weights, self.num_time_breaks-1, num_replicates))
        for i, x in enumerate(self.block_bootstrap(num_replicates, random_seed)):
            rates[:,:,i] = x
        return rates

