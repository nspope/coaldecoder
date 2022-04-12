# ------------ test fitting of demographic parameters -------------------#
# the example used is a big file, so we don't store it in this repo
# -----------------------------------------------------------------------#

library(coaldecoder)
reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

#---------------ancient samples------------------#

filename <- "asympt_iid_ancient_10x3.ts"
haps_per_pop <- 10
num_pops <- 3

#true parameters
time_breaks_true <- c(as.matrix(read.table(paste0(filename, ".t"))))
M_true <- array(c(as.matrix(read.table(paste0(filename, ".M")))), c(num_pops,num_pops,length(time_breaks_true)-1))
A_true <- array(c(as.matrix(read.table(paste0(filename, ".A")))), c(num_pops,num_pops,length(time_breaks_true)-1))

#transfer true parameters to finer time grid to compare
#observed/theoretical rates over time
time_breaks <- seq(0.0, max(time_breaks_true), 500)
fine_to_coarse <- as.numeric(cut(time_breaks, time_breaks_true))[-c(1)]
M <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
A <- array(0, c(num_pops,num_pops,length(time_breaks)-1))
for(i in 1:(length(time_breaks)-1))
{
  M[,,i] <- M_true[,,fine_to_coarse[i]]
  A[,,i] <- A_true[,,fine_to_coarse[i]]
}

pop_model <- PopulationTree$new("((A:30000,C:20000):10000,B:40000);", time_breaks=time_breaks)
pop_model$plot_population_tree()
pop_model$set_demographic_parameters(M)

# get observed rates
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- LETTERS[1:num_pops]
obs <- ObservedTrioRates(filename, sample_sets, time_breaks, bootstrap_blocks = 1000, bootstrap_replicates = 1000)

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

