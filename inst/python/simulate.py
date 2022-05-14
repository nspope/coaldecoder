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
import msprime
import tskit
import numpy as np
import pickle
import gzip

# TODO: add option not to concatenate but to write out separate ts

def msprime_simulate (sample_sizes, pop_times, pop_names, migr_mat, admix_mat, epoch_dur, length, recombination_rate, num_chromosomes, outfile, save, mutation_rate = 0.0, random_seed = None):

    assert save in ["tree_sequence", "msprime_inputs"]

    #TODO: checks

    num_pop = len(sample_sizes)
    num_epoch = len(epoch_dur)

    num_chromosomes = int(num_chromosomes)
    if random_seed is not None:
        random_seed = int(random_seed)

    assert migr_mat.shape[0] == num_pop and migr_mat.shape[1] == num_pop and migr_mat.shape[2] == num_epoch
    assert admix_mat.shape[0] == num_pop and admix_mat.shape[1] == num_pop and admix_mat.shape[2] == num_epoch
    assert len(pop_times) == num_pop
    assert length > 0
    assert recombination_rate >= 0.0
    assert num_chromosomes > 0
    assert len(pop_names) == num_pop

    # sample configurations
    populations = []
    for p in range(num_pop):
        populations += [msprime.PopulationConfiguration()]

    samples = []
    sample_names = []
    for p in range(num_pop):
        for s in range(int(sample_sizes[p])):
            samples += [(p, pop_times[p])]
            sample_names += [str(pop_names[p]) + str(s)]
    
    # demography
    demography = []
    t = 0
    for e in range(num_epoch):
        for i in range(num_pop):
            for j in range(num_pop):
                if i != j:
                    demography += [msprime.MassMigration(
                        time=t,
                        source=j,
                        dest=i,
                        proportion=admix_mat[i,j,e],
                    )]
        for i in range(num_pop):
            for j in range(num_pop):
                if i == j and np.isfinite(migr_mat[i,j,e]):
                    demography += [msprime.PopulationParametersChange(
                        time=t,
                        population=i, 
                        initial_size=migr_mat[i,j,e]/2,
                    )]
                elif i != j:
                    demography += [msprime.MigrationRateChange(
                        time=t,
                        matrix_index=(i, j), 
                        rate=migr_mat[i,j,e],
                    )]
        t += epoch_dur[e]

    if save == "msprime_inputs":
        msprime_inputs = {
            "demographic_events" : demography,
            "samples" : samples,
            "population_configurations" : populations,
            "recombination_rate" : recombination_rate,
            "mutation_rate" : mutation_rate,
            "length" : length,
            "num_replicates" : num_chromosomes,
            "random_seed" : random_seed,
        }
        pickle.dump(msprime_inputs, outfile)
    elif save == "tree_sequence":
        # simulation generator
        ts_gen = msprime.simulate(
          population_configurations=populations,
          samples=samples,
          demographic_events=demography,
          recombination_rate=recombination_rate,
          length=length,
          mutation_rate=mutation_rate,
          num_replicates=num_chromosomes,
          random_seed=random_seed,
        )
        
        # initialize
        ts = next(ts_gen)
        tables = ts.dump_tables().asdict()
        nodes = tables['nodes']
        edges = tables['edges']
        seq_offset = tables['sequence_length']
        node_offset = len(nodes['flags'])

        chromosome = 1
        if mutation_rate > 0.0:
            vcf_file = gzip.open(outfile + ".chr" + str(chromosome) + ".vcf.gz", "wt")
            ts.write_vcf(
                vcf_file, 
                ploidy=1, 
                individual_names=sample_names, 
                contig_id="chr" + str(chromosome)
            )
            vcf_file.close()
        
        # extract tables and offset values
        for ts in ts_gen:
            chromosome += 1
            if mutation_rate > 0.0:
                vcf_file = gzip.open(outfile + ".chr" + str(chromosome) + ".vcf.gz", "wt")
                ts.write_vcf(
                    vcf_file, 
                    ploidy=1, 
                    individual_names=sample_names, 
                    contig_id="chr" + str(chromosome)
                )
                vcf_file.close()

            tables = ts.dump_tables().asdict()
            nodes_flags = tables['nodes']['flags']
            nodes_time = tables['nodes']['time'][nodes_flags==0]
            nodes_offset = np.where(nodes_flags==0, node_offset - sum(nodes_flags==1), 0)
            edges_left = tables['edges']['left'] + seq_offset
            edges_right = tables['edges']['right'] + seq_offset
            edges_parent = tables['edges']['parent']
            edges_parent += nodes_offset[edges_parent]
            edges_child = tables['edges']['child']
            edges_child += nodes_offset[edges_child]
            
            # append to old table
            nodes['flags'] = np.append(nodes['flags'], [0] * len(nodes_time))
            nodes['time'] = np.append(nodes['time'], nodes_time)
            edges['left'] = np.append(edges['left'], edges_left)
            edges['right'] = np.append(edges['right'], edges_right)
            edges['parent'] = np.append(edges['parent'], edges_parent)
            edges['child'] = np.append(edges['child'], edges_child)
            
            # update offsets
            seq_offset += tables['sequence_length']
            node_offset += sum(nodes_flags==0)
        
        edge_sort = np.lexsort( (edges['left'], edges['child'], edges['parent'], nodes['time'][edges['parent']]) )
        tables = tskit.TableCollection(sequence_length=seq_offset)
        tables.nodes.set_columns(
            time=nodes['time'],
            flags=nodes['flags'].astype('uint32'),
        )
        tables.edges.set_columns(
            left=edges['left'][edge_sort],
            right=edges['right'][edge_sort],
            parent=edges['parent'][edge_sort],
            child=edges['child'][edge_sort],
        )
        ts = tables.tree_sequence()
        
        ts.dump(outfile + ".ts")

