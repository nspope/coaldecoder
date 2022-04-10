#-------------------------------------------------------------------------# 
# simulate data for three populations with population mergers,
# variable migration over time
#-------------------------------------------------------------------------#

import msprime
import tskit
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Simulate iid tree sequence with three populations")
parser.add_argument('--out', type=str, help='Path to output tree sequence')
parser.add_argument('--samples', type=int, help='Haploids per population')
parser.add_argument('--trees', type=int, help='Number of trees')
parser.add_argument('--seed', type=int, help='Random seed')
args = parser.parse_args()

# parameter arrays used by coaldecoder
M = np.zeros((3, 3, 6))
M[0,0,:] = 20000 * 2
M[0,1,0] = 1e-5; M[0,1,1] = 1e-4; M[0,1,2] = 1e-4; M[0,1,3] = 1e-5; M[0,1,4] = 1e-5; M[0,1,5:] = 0
M[0,2,0] = 0; M[0,2,1] = 1e-5; M[0,2,2] = 1e-6; M[0,2,3] = 1e-5; M[0,2,4:] = 0
M[1,1,:5] = 15000 * 2; M[1,1,5:] = M[0,0,5:]
M[1,0,:5] = 1e-5; M[1,0,5:] = 0
M[1,2,0] = 0; M[1,2,1:4] = 1e-5; M[1,2,4:] = 0
M[2,0,0] = 0; M[2,0,1:4] = 1e-5; M[2,0,4:] = 0
M[2,1,0] = 0; M[2,1,1:4] = 1e-5; M[2,1,4:] = 0
M[2,2,0] = np.Inf; M[2,2,1:4] = 10000 * 2; M[2,2,4:] = M[0,0,4:]

A = np.zeros((3, 3, 6))
A[0,0,:] = 1
A[1,1,:5] = 1; A[0,1,5:] = 1 #merger 1->0 at 40000
A[2,2,:4] = 1; A[0,2,4:] = 1 #merger 2->0 at 30000

t = np.array([0.0, 10000, 15000, 20000, 30000, 40000, 50000])

# write out flattened parameter arrays
M = M.flatten('F')
handle = open(args.out + ".M", 'w')
handle.write("# dim: 3 3 6\n")
for i in range(M.size): handle.write(str(M[i]) + "\n")

A = A.flatten('F')
handle = open(args.out + ".A", 'w')
handle.write("# dim: 3 3 6\n")
for i in range(A.size): handle.write(str(A[i]) + "\n")

handle = open(args.out + ".t", 'w')
handle.write("# dim: 7\n")
for i in range(t.size): handle.write(str(t[i]) + "\n")

# demography
populations = [
  msprime.PopulationConfiguration(initial_size=20000),
  msprime.PopulationConfiguration(initial_size=15000),
  msprime.PopulationConfiguration(initial_size=10000),
]

#super important -- whether or not to include acient population as "extant" from sampling point to present day
#yes: allow shit to migrate into it
#no: set migration rates to zero
#gives totally different answers
demography = [
  msprime.MigrationRateChange(time=0, rate=1e-5, matrix_index=(0,1)),
  msprime.MigrationRateChange(time=0, rate=0, matrix_index=(0,2)),
  msprime.MigrationRateChange(time=0, rate=1e-5, matrix_index=(1,0)),
  msprime.MigrationRateChange(time=0, rate=0, matrix_index=(1,2)),
  msprime.MigrationRateChange(time=0, rate=0, matrix_index=(2,0)),
  msprime.MigrationRateChange(time=0, rate=0, matrix_index=(2,1)),
  msprime.MigrationRateChange(time=10000, rate=1e-4, matrix_index=(0,1)),
  msprime.MigrationRateChange(time=10000, rate=1e-5, matrix_index=(0,2)),
  msprime.MigrationRateChange(time=10000, rate=1e-5, matrix_index=(1,0)),
  msprime.MigrationRateChange(time=10000, rate=1e-5, matrix_index=(1,2)),
  msprime.MigrationRateChange(time=10000, rate=1e-5, matrix_index=(2,0)),
  msprime.MigrationRateChange(time=10000, rate=1e-5, matrix_index=(2,1)),
  msprime.MigrationRateChange(time=15000, rate=1e-6, matrix_index=(0,2)),
  msprime.MigrationRateChange(time=20000, rate=1e-5),
  msprime.MassMigration(time=30000, source=2, dest=0),
  msprime.MigrationRateChange(time=30000, rate=0, matrix_index=(0,2)),
  msprime.MigrationRateChange(time=30000, rate=0, matrix_index=(1,2)),
  msprime.MigrationRateChange(time=30000, rate=0, matrix_index=(2,0)),
  msprime.MigrationRateChange(time=30000, rate=0, matrix_index=(2,1)),
  msprime.MassMigration(time=40000, source=1, dest=0),
  msprime.MigrationRateChange(time=40000, rate=0, matrix_index=(0,1)),
  msprime.MigrationRateChange(time=40000, rate=0, matrix_index=(1,0)),
]

# samples list
samples = []
for i in range(args.samples): samples += [(0, 0.0)]
for i in range(args.samples): samples += [(1, 0.0)]
for i in range(args.samples): samples += [(2, 10000.0)]

# simulation generator
ts_gen = msprime.simulate(
  population_configurations=populations,
  samples=samples,
  demographic_events=demography,
  Ne=20000,
  recombination_rate=0,
  length=1,
  num_replicates=args.trees,
  random_seed=args.seed
)

# initialize
ts = next(ts_gen)
tables = ts.dump_tables().asdict()
nodes = tables['nodes']
edges = tables['edges']
seq_offset = tables['sequence_length']
node_offset = len(nodes['flags'])

# extract tables and offset values
for ts in ts_gen:
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

print("Trees: " + str(ts.num_trees))
#print(ts.draw_text())

ts.dump(args.out)
