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

model = msprime.Demography.from_old_style(
    population_configurations=populations, 
    demographic_events=demography,
    ignore_sample_size=True,
)

ts = msprime.sim_ancestry(
  samples={i:3 for i in range(3)},
  demography=model,
  recombination_rate=1e-8,
  sequence_length=1e5,
  random_seed=1024,
  discrete_genome=True,
)

print("Trees: " + str(ts.num_trees))

ts.dump('example_3x3.ts')
