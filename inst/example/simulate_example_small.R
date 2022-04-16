#---------------------------------------------------------
# Simulate a very wacky-looking example for the purposes
# of illustration
#---------------------------------------------------------

library(coaldecoder)
#reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

num_pops <- 3
haps <- 10

time_breaks <- seq(0.0, 50000, length.out=1001)
epoch_start <- time_breaks[2:length(time_breaks)-1]
M <- array(0, c(num_pops,num_pops,length(epoch_start)))
M[1,1,] <- 100000 + 50000*cos(2*pi*epoch_start/15000)
M[1,2,] <- 1e-4 * exp(log(1e-6/1e-4)/40000 * epoch_start)
M[1,3,] <- 1e-5
M[2,1,] <- 1e-6 * exp(log(1e-4/1e-6)/40000 * epoch_start)
M[2,2,] <- 100000 + 50000*cos(2*pi*(epoch_start+15000)/15000)
M[2,3,] <- ifelse(epoch_start > 15000 & epoch_start < 25000, 1e-5, 1e-4)
M[3,1,] <- 1e-5
M[3,2,] <- ifelse(epoch_start > 15000 & epoch_start < 25000, 1e-5, 1e-6)
M[3,3,] <- 100000

pop_model <- PopulationTree$new("((A:30000,C:20000):10000,B:40000);", time_breaks=time_breaks)
pop_model$set_demographic_parameters(M)
pop_model$msprime_simulate(outfile="coaldecoder_example_small.ts", 
                           sample_sizes=c(haps,haps,haps), trees=1000, 
                           random_seed=1024, what="tree_sequence")

pop_assign <- rep(1:num_pops-1, each=haps)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- LETTERS[1:num_pops]

reticulate::source_python("../python/calculate_rates_multithreaded.py")
ptm_mt <- proc.time()
true_obs_rates <- ObservedTrioRatesMultithreaded("coaldecoder_example_small.ts", sample_sets, time_breaks, bootstrap_blocks=4, threads=2)
ptm_mt <- proc.time() - ptm_mt

reticulate::source_python("../python/calculate_rates.py")
ptm <- proc.time()
true_obs_rates_single <- ObservedTrioRates("coaldecoder_example_small.ts", sample_sets, time_breaks)
ptm <- proc.time() - ptm

all(dplyr::near(true_obs_rates_single$y, true_obs_rates$y))
