#---------------------------------------------------------
# Simulate a very wacky-looking example for the purposes
# of illustration
#---------------------------------------------------------

library(coaldecoder)
reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

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
pop_model$plot_demographic_parameters(time_scale=1000) + ggtitle("True demographic parameters") -> true_model_plot
pop_model$msprime_simulate(outfile="coaldecoder_example.ts", 
                           sample_sizes=c(haps,haps,haps), trees=50000, 
                           random_seed=1024, what="tree_sequence")

#pop_assign <- rep(1:num_pops-1, each=haps)
#sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
#names(sample_sets) <- LETTERS[1:num_pops]
#true_obs_rates <- ObservedTrioRates("coaldecoder_example.ts", sample_sets, time_breaks)
#
#true_rates <- true_obs_rates$rates()
#rownames(true_rates) <- true_obs_rates$emission_states()
#pop_model$plot_expected_coalescence_rates(log_transform=TRUE, time_scale=1000) + 
#  ggtitle("True vs simulated rates") -> true_obs_plot

#png("../inst/example/coaldecoder_example_true.png", height=7, width=14, units="in", res=300)
#cowplot::plot_grid(true_model, true_obs_plot)
#dev.off()

ggsave("../inst/example/coaldecoder_example_true.png", true_model_plot, height=7, width=7, units="in", dpi=300)

