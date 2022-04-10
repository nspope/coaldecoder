#-------------------------- "example_3x3.ts" --------------------------# 
# small example with three haps / three populations meant to test rate
# calculation with present-day samples
#----------------------------------------------------------------------#

import msprime

populations = [
  msprime.PopulationConfiguration(sample_size=3, initial_size=20000),
  msprime.PopulationConfiguration(sample_size=3, initial_size=10000),
  msprime.PopulationConfiguration(sample_size=3, initial_size=5000),
]

demography = [
  msprime.MigrationRateChange(time=0, rate=1e-5),
  msprime.MassMigration(time=30000, source=1, dest=0),
  msprime.MassMigration(time=30000, source=2, dest=0),
]

ts = msprime.simulate(
  population_configurations=populations,
  demographic_events=demography,
  Ne=20000,
  recombination_rate=1e-8,
  length=1e5,
  random_seed=1024
)

print("Trees: " + str(ts.num_trees))

ts.dump('example_3x3.ts')
