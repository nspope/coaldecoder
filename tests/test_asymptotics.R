# ------------ test asymptotic correctness of rates ---------------------#
# e.g. concordence between observed/expected rates for lots of iid trees
# the example used is a big file, so we don't store it in this repo
# it needs to be generated via 'python simulate_3pop_iid.py' and
# 'python simulate_3pop_iid_ancient.py'
#
# python3 simulate_3pop_iid_ancient.py --samples 10 --trees 50000 --out asympt_iid_ancient_10x3.ts --seed 1024
# python3 simulate_3pop_iid.py --samples 10 --trees 50000 --out asympt_iid_10x3.ts --seed 1024
# -----------------------------------------------------------------------#

library(coaldecoder)
reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

# -------------------- all contemporary samples -------------------------#

{
filename <- "asympt_iid_10x3.ts"
haps_per_pop <- 10
num_pops <- 3

#true parameters
time_breaks_true <- c(as.matrix(read.table(paste0(filename, ".t"))))
M_true <- array(c(as.matrix(read.table(paste0(filename, ".M")))), c(num_pops,num_pops,length(time_breaks_true)-1))
A_true <- array(c(as.matrix(read.table(paste0(filename, ".A")))), c(num_pops,num_pops,length(time_breaks_true)-1))

#transfer true parameters to finer time grid to compare
#observed/theoretical rates over time
time_breaks <- seq(0.0, max(time_breaks_true), 2500)
fine_to_coarse <- as.numeric(cut(time_breaks, time_breaks_true))[-c(1)]
M <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
A <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
for(i in 1:(length(time_breaks)-1))
{
  M[,,i] <- M_true[,,fine_to_coarse[i]]
  A[,,i] <- A_true[,,fine_to_coarse[i]]
}

#TODO: separate set of tests for these members
pop_model <- PopulationTree$new("((A:30000,C:30000):10000,B:40000);", time_breaks=time_breaks)
pop_model$plot_population_tree()
pop_model$set_demographic_parameters(M)
pop_model$demographic_parameters()
pop_model$admixture_coefficients()
pop_model$time_breaks()
pop_model$epoch_durations()
pop_model$population_sampling_times()
pop_model$plot_expected_coalescence_rates()

# get observed rates
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- LETTERS[1:num_pops]
obs <- ObservedTrioRates(filename, sample_sets, time_breaks)

# get theoretical rates
exp_rates <- pop_model$expected_coalescence_rates()

# check against original parameter arrays
deco <- CoalescentDecoder$new(num_pops, diff(time_breaks), TRUE)
stopifnot(all(dplyr::near(deco$expected_rates(deco$initial_state_vectors(), M, A)$y, exp_rates)))

# format observed rates
obs_rates <- obs$rates()
rownames(obs_rates) <- obs$emission_states()
colnames(obs_rates) <- 1:ncol(obs_rates)
stopifnot(rownames(exp_rates) == rownames(obs_rates))

# plot rates
pop_model$plot_expected_coalescence_rates(observed=obs_rates, log_transform=TRUE, time_scale=1000) + 
  ggtitle("True vs simulated rates")
}

# -------------------- one population has ancient samples -------------------------#

{
filename <- "asympt_iid_ancient_10x3.ts"
haps_per_pop <- 10
num_pops <- 3

#true parameters
time_breaks_true <- c(as.matrix(read.table(paste0(filename, ".t"))))
M_true <- array(c(as.matrix(read.table(paste0(filename, ".M")))), c(num_pops,num_pops,length(time_breaks_true)-1))
A_true <- array(c(as.matrix(read.table(paste0(filename, ".A")))), c(num_pops,num_pops,length(time_breaks_true)-1))

#transfer true parameters to finer time grid to compare
#observed/theoretical rates over time
time_breaks <- seq(0.0, max(time_breaks_true), 2500)
fine_to_coarse <- as.numeric(cut(time_breaks, time_breaks_true))[-c(1)]
M <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
A <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
for(i in 1:(length(time_breaks)-1))
{
  M[,,i] <- M_true[,,fine_to_coarse[i]]
  A[,,i] <- A_true[,,fine_to_coarse[i]]
}

#TODO: separate set of tests for these members
pop_model <- PopulationTree$new("((A:30000,C:20000):10000,B:40000);", time_breaks=time_breaks)
pop_model$plot_population_tree()
pop_model$set_demographic_parameters(M)
pop_model$demographic_parameters()
pop_model$admixture_coefficients()
pop_model$time_breaks()
pop_model$epoch_durations()
pop_model$population_sampling_times()
pop_model$plot_expected_coalescence_rates()

# get observed rates
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- LETTERS[1:num_pops]
obs <- ObservedTrioRates(filename, sample_sets, time_breaks)

# get theoretical rates
exp_rates <- pop_model$expected_coalescence_rates()

# check against original parameter arrays
deco <- CoalescentDecoder$new(num_pops, diff(time_breaks), TRUE)
stopifnot(all(dplyr::near(deco$expected_rates(deco$initial_state_vectors(), M, A)$y, exp_rates)))

# format observed rates
obs_rates <- obs$rates()
rownames(obs_rates) <- obs$emission_states()
colnames(obs_rates) <- 1:ncol(obs_rates)
stopifnot(rownames(exp_rates) == rownames(obs_rates))

# plot rates
pop_model$plot_expected_coalescence_rates(observed=obs_rates, log_transform=TRUE, time_scale=1000) + 
  ggtitle("True vs simulated rates")

#####-----------######

pop_model$plot_expected_coalescence_rates(observed=obs_rates, log_transform=TRUE, time_scale=1000) + 
  ggtitle("True vs simulated rates") -> rates_plot
pop_model$plot_demographic_parameters(time_scale=1000) + 
  ggtitle("Demographic model") -> model_plot

pdf("example_ancient_samples.pdf", height=7, width=14)
cowplot::plot_grid(model_plot, rates_plot)
dev.off()

}

#-----------------------check simulation validity-------------------#
#we use the same parameters as in simulate_3pop_iid_ancient.py
{
filename <- "asympt_iid_ancient_10x3.ts"
haps_per_pop <- 10
num_pops <- 3

#true parameters
time_breaks_true <- c(as.matrix(read.table(paste0(filename, ".t"))))
M_true <- array(c(as.matrix(read.table(paste0(filename, ".M")))), c(num_pops,num_pops,length(time_breaks_true)-1))
A_true <- array(c(as.matrix(read.table(paste0(filename, ".A")))), c(num_pops,num_pops,length(time_breaks_true)-1))

#transfer true parameters to finer time grid to compare
#observed/theoretical rates over time
time_breaks <- seq(0.0, max(time_breaks_true), 2500)
fine_to_coarse <- as.numeric(cut(time_breaks, time_breaks_true))[-c(1)]
M <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
A <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
for(i in 1:(length(time_breaks)-1))
{
  M[,,i] <- M_true[,,fine_to_coarse[i]]
  A[,,i] <- A_true[,,fine_to_coarse[i]]
}

pop_model <- PopulationTree$new("((A:30000,C:20000):10000,B:40000);", time_breaks=time_breaks)
pop_model$set_demographic_parameters(M)

pop_model$msprime_simulate(outfile=paste0(filename, ".check"), sample_sizes=c(10,10,10), trees=50000, random_seed=1024, what="tree_sequence")

#check observed rates
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- LETTERS[1:num_pops]
obs2 <- ObservedTrioRates(paste0(filename, ".check"), sample_sets, time_breaks)

obs_rates <- obs2$rates()
rownames(obs_rates) <- obs2$emission_states()
pop_model$plot_expected_coalescence_rates(observed=obs_rates, log_transform=TRUE, time_scale=1000) + 
  ggtitle("True vs simulated rates")

}
