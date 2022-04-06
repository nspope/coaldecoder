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
import msprime
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

    if isinstance(ts_paths, str):
        ts_paths = [ts_paths]

    sample_sets_list = []
    for p in sample_sets.keys():
        sample_sets_list.append([int(i) for i in sample_sets[p]])
    sample_sets = sample_sets_list

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
                num_replicates=bootstrap_replicates
            )
            for distr_boot, i in zip(generator, range(bootstrap_replicates)):
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

def _count_all_trio_coalescence_events(node, tree, sample_sets):
    """
    TODO: this also counts first coalescence events, see above.
    Version that makes it into tskit should be *only* full trios.

    Count the number of trios that coalesce in node,
    within and between the sets of samples in ``sample_sets``. In other
    words, count topologies of the form ``((A,B),C)node`` where ``A,B,C``
    are labels and `node` is the node ID.  The count of pairs with members
    that belong to sets :math:`a` and :math:`b` with outgroup :math:`c` is:

    .. math:

        TODO

    where :math:`C_i(a)` is the number of samples from set :math:`a`
    descended from child :math:`i` of the node, and :math:`O(c)` is the
    number of samples from set :math:`c` that are *not* descended from the
    node.  The values in the output are ordered canonically by pair then
    outgroup; e.g. if ``len(sample_sets) == 2`` then the values would
    correspond to counts of pairs with set labels,
    ``[((0,0),0), ((0,0),1), ((0,1),0), ((0,1),1), ...]``.
    """
    samples = list(tree.samples(node))
    outg_counts = [len(s) - len(np.intersect1d(samples, s)) for s in sample_sets]
    pair_counts = CoalescenceTimeDistribution._count_pair_coalescence_events(
        node, tree, sample_sets
    )
    trio_counts = []
    # first coalescence events
    for i in pair_counts:
        for j in outg_counts:
            trio_counts.append(i * j)
    # second coalescence events
    children = tree.children(node)
    child_pair_counts = np.zeros((len(children), len(pair_counts)), "int32")
    child_sample_counts = np.zeros((len(children), len(outg_counts)), "int32")
    for i,child in enumerate(children):
        samples = list(tree.samples(child))
        child_sample_counts[i,:] = [
            len(np.intersect1d(samples, s)) for s in sample_sets
        ]
        child_pair_counts[i,:] = CoalescenceTimeDistribution._count_pair_coalescence_events(
            child, tree, sample_sets
        )
    #for i in possible pairs:
    #  for j in possible samples:
    #       sum = 0
    #       iterate over children duos u,v
    #           sum += number of i-pair in u * number of j-sample in v +
    #                  number of i-pair in v * number of j-sample in u
    for i in range(child_pair_counts.shape[1]):
        for j in range(child_sample_counts.shape[1]):
            count = 0
            for u, v in itertools.combinations(
                range(len(children)), 2
            ):
                count += (child_pair_counts[u,i] * child_sample_counts[v,j] +
                    child_pair_counts[v,i] * child_sample_counts[u,j])
            trio_counts.append(count)
    return np.array(trio_counts, "int32")

def trio_coalescence_rates(
    ts_paths,
    sample_sets,
    time_breaks,
    bootstrap_replicates = 0,
    bootstrap_blocks = 1,
    random_seed = None,
    ):

    if isinstance(ts_paths, str):
        ts_paths = [ts_paths]

    sample_sets_list = []
    for p in sample_sets.keys():
        sample_sets_list.append([int(i) for i in sample_sets[p]])
    sample_sets = sample_sets_list

    time_breaks = np.append(np.array(time_breaks, dtype="float"), 0)
    time_breaks = np.sort(np.unique(time_breaks))

    bootstrap_replicates = int(bootstrap_replicates)
    bootstrap_blocks = int(bootstrap_blocks)
    if random_seed is not None:
        random_seed = int(random_seed)

    out = []
    for path in ts_paths:
        # read tree sequence
        ts = tskit.load(path)
        distr = CoalescenceTimeDistribution(
            ts,
            weight_func=_count_all_trio_coalescence_events, 
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
                num_replicates=bootstrap_replicates
            )
            for distr_boot, i in zip(generator, range(bootstrap_replicates)):
                n_tmp = np.nan_to_num(distr_boot.num_uncoalesced(np.array([0.0])))
                n_boot[:,i] = n_tmp[:,0,0]
                y_tmp = np.nan_to_num(distr_boot.num_coalesced(time_breaks))
                y_boot[:,:,i] = y_tmp[:,1:,0] - y_tmp[:,:-1,0]
        out += [{
            "file" : path, "y" : y, "n" : n, 
            "y_boot" : y_boot, "n_boot" : n_boot, 
            "size" : ts.sequence_length, "trees" : ts.num_trees, 
        }]
    return out

#def test_TrioCounter(n=30, seed=1):
#    msprime.simulate(, random_seed=1)
#    #itertools get unique trios, make map
#    for a,b in itertools.combinations_with_replacement()
#        trio_map[str("
#    #loop over trios
#    #count number of labellings for each node


#----------------- MOCKUP FOR RATE CALCULATOR
# this is the basis for https://github.com/tskit-dev/tskit/pull/2119 
# DELETME when that is merged into dev branch

class CoalescenceTimeTable:
    """
    Container for sorted coalescence times, weights, and block assignments.
    """

    def __init__(self, time, block, weights):
        assert time.shape[0] == weights.shape[0] == block.shape[0]
        assert time.ndim == block.ndim == 1
        assert weights.ndim == 2
        self.num_weights = weights.shape[1]
        # remove empty records
        not_empty = np.sum(weights, 1) > 0
        self.num_records = sum(not_empty)
        self.time = time[not_empty]
        self.block = block[not_empty]
        self.weights = weights[not_empty, :]
        # add left boundary at time 0
        self.num_records += 1
        self.time = np.pad(self.time, (0, 1))
        self.block = np.pad(self.block, (0, 1))
        self.weights = np.pad(self.weights, ((0, 1), (0, 0)))
        # sort by node time
        time_order = np.argsort(self.time)
        self.time = self.time[time_order]
        self.block = self.block[time_order]
        self.weights = self.weights[time_order, :]
        # calculate quantiles
        self.num_blocks = 1 + np.max(self.block) if self.num_records > 0 else 0
        self.block_multiplier = np.ones(self.num_blocks)
        self.cum_weights = np.cumsum(self.weights, 0)
        self.quantile = np.empty((self.num_records, self.num_weights))
        self.quantile[:] = np.nan
        for i in range(self.num_weights):
            if self.cum_weights[-1, i] > 0:
                self.quantile[:, i] = self.cum_weights[:, i] / self.cum_weights[-1, i]

    def resample_blocks(self, block_multiplier):
        assert block_multiplier.shape[0] == self.num_blocks
        assert np.sum(block_multiplier) == self.num_blocks
        self.block_multiplier = block_multiplier
        for i in range(self.num_weights):
            self.cum_weights[:, i] = np.cumsum(
                self.weights[:, i] * self.block_multiplier[self.block], 0
            )
            if self.cum_weights[-1, i] > 0:
                self.quantile[:, i] = self.cum_weights[:, i] / self.cum_weights[-1, i]


class CoalescenceTimeDistribution:
    """
    Class to precompute a table of sorted/weighted node times, from which to calculate
    the empirical distribution function and estimate coalescence rates in time windows.
    """

    @staticmethod
    def _count_coalescence_events(node, tree, sample_sets):
        return np.array([1], dtype="int32")

    @staticmethod
    def _count_pair_coalescence_events(node, tree, sample_sets):
        """
        Count the number of pairs that coalesce in node, within and between the
        sets of samples in ``sample_sets``. The count of pairs with members that
        belong to sets :math:`a` and :math:`b` is:

        .. math:

            \\sum_{i \\neq j} (C_i(a) C_j(b) + C_i(b) C_j(a))/(1 - \\mathbb{I}[a = b])

        where :math:`C_i(a)` is the number of samples from set :math:`a`
        descended from child :math:`i`.  The values in the output are ordered
        canonically; e.g. if ``len(sample_sets) == 2`` then the values would
        correspond to counts of pairs with set labels ``[(0,0), (0,1), (1,1)]``.
        """

        # TODO needs to be optimized, use np.intersect1d
        children = tree.children(node)
        samples_per_child = [set(list(tree.samples(c))) for c in children]
        sample_counts = np.zeros((len(sample_sets), len(children)), "int32")
        for i, s1 in enumerate(samples_per_child):
            for a, s2 in enumerate([set(s) for s in sample_sets]):
                sample_counts[a, i] = len(s1 & s2)

        pair_counts = []
        for a, b in itertools.combinations_with_replacement(
            range(sample_counts.shape[0]), 2
        ):
            count = 0
            for i, j in itertools.combinations(range(sample_counts.shape[1]), 2):
                count += (
                    sample_counts[a, i] * sample_counts[b, j]
                    + sample_counts[a, j] * sample_counts[b, i]
                ) / (1 + int(a == b))
            pair_counts.append(count)

        return np.array(pair_counts, "int32")

    @staticmethod
    def _count_trio_first_coalescence_events(node, tree, sample_sets):
        """
        Count the number of pairs that coalesce in node with an outgroup,
        within and between the sets of samples in ``sample_sets``. In other
        words, count topologies of the form ``((A,B):node,C)`` where ``A,B,C``
        are labels and `node` is the node ID.  The count of pairs with members
        that belong to sets :math:`a` and :math:`b` with outgroup :math:`c` is:

        .. math:

            \\sum_{i \\neq j} (C_i(a) C_j(b) + C_i(b) C_j(a)) \\times
            O(c) / (1 - \\mathbb{I}[a = b])

        where :math:`C_i(a)` is the number of samples from set :math:`a`
        descended from child :math:`i` of the node, and :math:`O(c)` is the
        number of samples from set :math:`c` that are *not* descended from the
        node.  The values in the output are ordered canonically by pair then
        outgroup; e.g. if ``len(sample_sets) == 2`` then the values would
        correspond to counts of pairs with set labels,
        ``[((0,0),0), ((0,0),1), ((0,1),0), ((0,1),1), ...]``.
        """
        samples = list(tree.samples(node))
        outg_counts = [len(s) - len(np.intersect1d(samples, s)) for s in sample_sets]
        pair_counts = CoalescenceTimeDistribution._count_pair_coalescence_events(
            node, tree, sample_sets
        )
        trio_counts = []
        for i in pair_counts:
            for j in outg_counts:
                trio_counts.append(i * j)
        return np.array(trio_counts, "int32")

    def _update_weights_by_edge_diff(self, tree, edge_diff, running_weights):
        """
        Update ``running_weights`` to reflect ``tree`` using edge differences
        ``edge_diff`` with the previous tree.
        """

        assert edge_diff.interval == tree.interval

        # nodes that have been removed from tree
        removed = {i.child for i in edge_diff.edges_out if tree.is_isolated(i.child)}
        # TODO: What if sample is removed from tree? In that case should all
        # nodes be updated?

        # nodes where descendant subtree has been altered
        modified = {i.parent for i in edge_diff.edges_in}
        for i in copy.deepcopy(modified):
            while not tree.parent(i) is tskit.NULL and not tree.parent(i) in modified:
                i = tree.parent(i)
                modified.add(i)

        # recalculate weights for current tree
        for i in removed:
            running_weights[i, :] = 0
        for i in modified:
            running_weights[i, :] = self.weight_func(i, tree, self.sample_sets)
        self.weight_func_evals += len(modified)

    def _build_ecdf_table_for_window(
        self, left, right, tree, edge_diffs, running_weights
    ):
        """
        Construct ECDF table for genomic interval [left, right]. Update ``tree``,
        ``edge_diffs``, and ``running_weights`` for input for next window. Trees are
        counted as belonging to any interval with which they overlap, and thus
        can be used in several intervals. Thus, the concatenation of ECDF
        tables across multiple intervals is not the same as the ECDF table
        for the union of those intervals. Trees within intervals are chunked
        into roughly equal-sized blocks for bootstrapping.
        """

        assert tree.interval.left <= left and right > left

        # assign trees in window to equal-sized blocks with unique id
        _tree = tree.copy()
        if right >= _tree.tree_sequence.sequence_length:
            _tree.last()
        else:
            # _tree.seek(right) won't work if `right` is recomb breakpoint
            while _tree.interval.right < right:
                _tree.next()
        tree_idx = np.arange(tree.index, _tree.index + 1) - tree.index
        tree_offset = tree.index
        num_blocks = min(self.num_blocks, len(tree_idx))
        tree_blocks = np.floor_divide(num_blocks * tree_idx, len(tree_idx))

        # calculate span weights
        _tree = tree.copy()
        tree_span = [
            min(_tree.interval.right, right) - 
            max(_tree.interval.left, left)
        ]
        while _tree.index < tree_offset + tree_idx[-1]:
            _tree.next()
            tree_span.append(
                min(_tree.interval.right, right) - 
                max(_tree.interval.left, left)
            )
        tree_span = np.array(tree_span) / sum(tree_span)

        # storage if using single window, block for entire tree sequence
        buffer_size = self.buffer_size
        table_size = buffer_size
        time = np.zeros(table_size)
        block = np.zeros(table_size, "int32")
        weights = np.zeros((table_size, self.num_weights))

        # assemble table of coalescence times in window
        indices = np.zeros(tree.tree_sequence.num_nodes, "int32") - 1
        last_block = np.zeros(tree.tree_sequence.num_nodes, "int32") - 1
        num_record = 0
        while tree != tskit.NULL:
            if tree.interval.right > left:
                current_block = tree_blocks[tree.index - tree_offset]
                if self.span_normalise:
                    span_weight = tree_span[tree.index - tree_offset]
                else:
                    span_weight = 1.0
                nodes_in_tree = np.array(
                    [i for i in tree.nodes() if tree.is_internal(i)]
                )
                # TODO nodes_in_tree should take into account isolated samples,
                # and behaviour below should be adjusted
                nodes_to_add = nodes_in_tree[
                    np.where(last_block[nodes_in_tree] != current_block)
                ]
                if len(nodes_to_add) > 0:
                    idx = np.arange(num_record, num_record + len(nodes_to_add))
                    last_block[nodes_to_add] = current_block
                    indices[nodes_to_add] = idx
                    if table_size < num_record + len(nodes_to_add):
                        table_size += buffer_size
                        time = np.pad(time, (0, buffer_size))
                        block = np.pad(block, (0, buffer_size))
                        weights = np.pad(weights, ((0, buffer_size), (0, 0)))
                    time[idx] = [tree.time(i) for i in nodes_to_add]
                    block[idx] = current_block
                    num_record += len(nodes_to_add)
                weights[indices[nodes_in_tree], :] += span_weight * running_weights[nodes_in_tree, :]

            if tree.interval.right < right:
                # if current tree does not cross window boundary, move to next
                tree.next()
                self._update_weights_by_edge_diff(
                    tree, next(edge_diffs), running_weights
                )
            else:
                # use current tree as initial tree for next window
                break

        return CoalescenceTimeTable(time, block, weights)

    def _generate_ecdf_tables(self, ts, window_breaks):
        """
        Return generator for ECDF tables across genomic windows defined by
        ``window_breaks``.

        ..note:: This could be used in methods in place of loops over
            pre-assembled tables.
        """

        tree = ts.first()
        edge_diffs = ts.edge_diffs()
        running_weights = np.zeros((ts.num_nodes, self.num_weights))
        self._update_weights_by_edge_diff(tree, next(edge_diffs), running_weights)
        for left, right in zip(window_breaks[:-1], window_breaks[1:]):
            yield self._build_ecdf_table_for_window(
                left, right, tree, edge_diffs, running_weights
            )

    def __init__(
        self,
        ts,
        sample_sets=None,
        weight_func=None,
        window_breaks=None,
        blocks_per_window=None,
        span_normalise=True,
    ):

        assert isinstance(ts, tskit.trees.TreeSequence)

        if sample_sets is None:
            sample_sets = [list(ts.samples())]
        assert all([isinstance(i, list) for i in sample_sets])
        assert all([i in ts.samples() for j in sample_sets for i in j])
        self.sample_sets = sample_sets

        if weight_func is None or weight_func == "coalescence_events":
            self.weight_func = self._count_coalescence_events
        elif weight_func == "pair_coalescence_events":
            self.weight_func = self._count_pair_coalescence_events
        elif weight_func == "trio_first_coalescence_events":
            self.weight_func = self._count_trio_first_coalescence_events
        else:
            assert callable(weight_func)
            self.weight_func = weight_func
        _weight_func_eval = self.weight_func(0, ts.first(), self.sample_sets)
        assert isinstance(_weight_func_eval, np.ndarray)
        assert _weight_func_eval.ndim == 1
        self.num_weights = len(_weight_func_eval)

        if window_breaks is None:
            window_breaks = np.array([0.0, ts.sequence_length])
        assert isinstance(window_breaks, np.ndarray)
        assert window_breaks.ndim == 1
        assert np.min(window_breaks) >= 0.0
        assert np.max(window_breaks) <= ts.sequence_length
        window_breaks = np.sort(np.unique(window_breaks))
        self.windows = [
            tskit.trees.Interval(left, right)
            for left, right in zip(window_breaks[:-1], window_breaks[1:])
        ]
        self.num_windows = len(self.windows)

        if blocks_per_window is None:
            blocks_per_window = 1
        assert isinstance(blocks_per_window, int)
        assert blocks_per_window > 0
        self.num_blocks = blocks_per_window

        assert isinstance(span_normalise, bool)
        self.span_normalise = span_normalise

        self.buffer_size = ts.num_nodes
        self.weight_func_evals = 0
        self.tables = [table for table in self._generate_ecdf_tables(ts, window_breaks)]

    # TODO
    #
    # def __str__(self):
    #    return self.useful_text_summary()
    #
    # def __repr_html__(self):
    #    return self.useful_html_summary()

    def copy(self):
        return copy.deepcopy(self)

    def ecdf(self, times):
        """
        Returns the empirical distribution function evaluated at the time
        points in ``times``.

        The output array has shape ``(self.num_weights, len(times),
        self.num_windows)``.
        """

        assert isinstance(times, np.ndarray)
        assert times.ndim == 1

        values = np.empty((self.num_weights, len(times), self.num_windows))
        values[:] = np.nan
        for k, table in enumerate(self.tables):
            indices = np.searchsorted(table.time, times, side="right") - 1
            assert all([0 <= i < table.num_records for i in indices])
            values[:, :, k] = table.quantile[indices, :].T
        return values

    # TODO
    #
    # def quantile(self, times):
    #   """
    #   Return interpolated quantiles of coalescence times, using the same
    #   approach as numpy.quantile(..., method="linear")
    #   """

    def num_coalesced(self, times):
        """
        Returns number of coalescence events that have occured by the time
        points in ``times``.

        The output array has shape ``(self.num_weights, len(times),
        self.num_windows)``.
        """

        assert isinstance(times, np.ndarray)
        assert times.ndim == 1

        values = self.ecdf(times)
        for k, table in enumerate(self.tables):
            weight_totals = table.cum_weights[-1, :].reshape(values.shape[0], 1)
            values[:, :, k] *= np.tile(weight_totals, (1, values.shape[1]))
        return values

    def num_uncoalesced(self, times):
        """
        Returns the number of coalescence events remaining by the time points
        in ``times``.

        The output array has shape ``(self.num_weights, len(times),
        self.num_windows)``.
        """

        values = 1.0 - self.ecdf(times)
        for k, table in enumerate(self.tables):
            weight_totals = table.cum_weights[-1, :].reshape(values.shape[0], 1)
            values[:, :, k] *= np.tile(weight_totals, (1, values.shape[1]))
        return values

    def mean(self, since=0.0):
        """
        Returns the average time between ``since`` and the coalescence events
        that occurred after ``since``.

        Note that ``1/self.mean(left)`` is an estimate of the coalescence rate
        over the interval (left, infinity).

        The output array has shape ``(self.num_weights, self.num_windows)``.

        ..note:: Check for overflow in ``np.average``.
        """

        assert isinstance(since, float) and since >= 0.0

        values = np.empty((self.num_weights, self.num_windows))
        values[:] = np.nan
        for k, table in enumerate(self.tables):
            index = np.searchsorted(table.time, since, side="right")
            if index == table.num_records:
                values[:, k] = np.nan
            else:
                for i in range(self.num_weights):
                    if table.cum_weights[-1, i] > 0:
                        multiplier = table.block_multiplier[table.block[index:]]
                        values[i, k] = np.average(
                            table.time[index:] - since,
                            weights=table.weights[index:, i] * multiplier,
                        )
        return values

    def coalescence_probability_in_intervals(self, time_breaks):
        """
        Returns the proportion of coalescence events occurring in the time
        intervals defined by ``time_breaks``, out of events that have not
        yet occurred by the intervals' left boundaries.

        The output array has shape ``(self.num_weights, len(time_breaks)-1,
        self.num_windows)``.
        """

        assert isinstance(time_breaks, np.ndarray)

        time_breaks = np.sort(np.unique(time_breaks))
        num_coalesced = self.num_coalesced(time_breaks)
        num_uncoalesced = self.num_uncoalesced(time_breaks)
        numer = num_coalesced[:, 1:, :] - num_coalesced[:, :-1, :]
        denom = num_uncoalesced[:, :-1, :]
        return numer / np.where(np.isclose(denom, 0.0), np.nan, denom)

    def coalescence_rate_in_intervals(self, time_breaks):
        """
        Returns the interval-censored Kaplan-Meier estimate of the hazard rate for
        coalesence events within the time intervals defined by ``time_breaks``. The
        estimator is,

        .. math::
            \\hat{c}_{l,r} = \\begin{cases}
              \\log(1 - x_{l,r}/k_{l})/(l - r) & \\mathrm{if~} x_{l,r} < k_{l} \\\\
              \\hat{c}_{l,r} = k_{l} / t_{l,r} & \\mathrm{if~} x_{l,r} = k_{l}
            \\end{cases}

        and is undefined where :math:`k_{l} = 0`. Here, :math:`x_{l,r}` is the
        number of events occuring in time interval :math:`(l, r]`,
        :math:`k_{l}` is the number of events remaining at time :math:`l`, and
        :math:`t_{l,r}` is the sum of event times occurring in the interval
        :math:`(l, r]`.

        The output array has shape ``(self.num_weights, len(time_breaks)-1,
        self.num_windows)``.
        """

        assert isinstance(time_breaks, np.ndarray)

        time_breaks = np.sort(np.unique(time_breaks))
        phi = self.coalescence_probability_in_intervals(time_breaks)
        duration = np.reshape(time_breaks[1:] - time_breaks[:-1], (1, phi.shape[1], 1))
        numer = -np.log(1.0 - np.where(np.isclose(phi, 1.0), np.nan, phi))
        denom = np.tile(duration, (self.num_weights, 1, self.num_windows))
        for i, j, k in np.argwhere(np.isclose(phi, 1.0)):
            numer[i, j, k] = 1.0
            denom[i, j, k] = self.mean(time_breaks[j])[i, k]
        return numer / denom

    def block_bootstrap(self, num_replicates=1, random_seed=None):
        """
        Return a generator that produces ``num_replicates`` copies of the
        object where blocks within genomic windows are randomly resampled.

        ..note:: Copying could be expensive.
        """

        rng = np.random.default_rng(random_seed)
        for _i in range(num_replicates):
            replicate = self.copy()
            for table in replicate.tables:
                block_multiplier = rng.multinomial(
                    table.num_blocks, [1.0 / table.num_blocks] * table.num_blocks
                )
                table.resample_blocks(block_multiplier)
            yield replicate
