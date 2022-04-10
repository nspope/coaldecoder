# ------------ test asymptotic correctness of rates ---------------------#
# e.g. concordence between observed/expected rates for lots of iid trees
# the example used is a big file, so we don't store it in this repo
# it needs to be generated via 'python asympt_iid_3x3.py'
# -----------------------------------------------------------------------#

library(coaldecoder)
reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))
reticulate::source_python("../inst/python/calculate_rates.py")

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

# get observed rates
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
obs <- ObservedTrioRates(filename, sample_sets, time_breaks)

# get theoretical rates
deco <- CoalescentDecoder$new(num_pops, diff(time_breaks), TRUE)
ll <- deco$loglikelihood(deco$initial_state_vectors(), M, A)
exp_rates <- ll$y_hat

# plot rates
obs_rates <- obs$rates()
rownames(obs_rates) <- obs$emission_states()
colnames(obs_rates) <- 1:ncol(obs_rates)
plot_coalescence_rates(obs_rates, exp_rates)

# plot parameters
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

# get observed rates
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
obs <- ObservedTrioRates(filename, sample_sets, time_breaks)

# get theoretical rates
deco <- CoalescentDecoder$new(num_pops, diff(time_breaks), TRUE)
ll <- deco$loglikelihood(deco$initial_state_vectors(), M, A)
exp_rates <- ll$y_hat

# plot rates
obs_rates <- obs$rates()
rownames(obs_rates) <- obs$emission_states()
colnames(obs_rates) <- 1:ncol(obs_rates)
rates_plot <- plot_coalescence_rates(obs_rates, diff(time_breaks), exp_rates, time_scale=1000, log_transform=TRUE) + ggtitle("True vs simulated rates")

# plot parameters
model_plot <- plot_demographic_model(M, diff(time_breaks), time_scale=1000) + ggtitle("Demographic model")

pdf("example_ancient_samples.pdf", height=7, width=14)
cowplot::plot_grid(model_plot, rates_plot)
dev.off()

}

#---------------DELETEME
# API thoughts
#
# DemographicModel(time_breaks=time_breaks,
#                  populations = list(
#                     list(population_id, sampling_time, merges_into=c(destination, epoch),
#                     list(population_id, sampling_time, merges_into=c(destination, epoch)
#                 )
# )
# DemographicModel$get_demographic_parameters()
# DemographicModel$set_effective_population_size(index, which_epochs, value)
# DemographicModel$set_migration_rate(matrix_index, which_epochs, value)
# DemographicModel$simulate("example.ts", samples = c(10, 10, 10), recombination_rate = "iid", size = 1000)
