# ---------------------------------------------------------------------------- #
# 11 May 22                                                                    #
# Validation simulations for demographic inference from trio coalescence rates #
# Two populations, all contemporary samples, recombinant sequences             #
# ---------------------------------------------------------------------------- #

#devtools::install_github("nspope/coaldecoder")
library(coaldecoder)
reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

# ------------- settings ------------ #
random_seed <- 1025                              #reproducibility will be system-specific, sadly
haploids_per_population <- 10                    #number of haploid samples per population
max_time <- 25000                                #time until which to estimate demographic parameters
number_of_epochs <- 25                           #number of epochs across which to estimate demographic parameters
number_of_chromosomes <- 30                      #number of chromosomes to use to estimate demographic parameters
holdout_chromosomes <- 10                        #number of chromosomes to use as holdout for selecting smoothing penalty
chromosome_length <- 10e6                        #chromosome length in base pairs
recombination_rate <- 1e-8                       #recombination rate; if zero there is a single tree per chromosome
mutation_rate <- 1e-8                            #per base mutation rate; if zero no VCF are produced
output_dir <- paste0('simulation_', random_seed) #directory to dump files into
output_prefix <- paste0(output_dir, '/sim')      #prefix to use for output files
smoothing_penalty <- seq(0, 30, 1)               #fit models across this grid of penalties

#generate IID genealogies instead of recombinant tree sequences 
#(much faster, less realistic, useful for checking identifiability)
iid_genealogies <- FALSE
if (iid_genealogies){
  number_of_chromosomes <- 20000
  holdout_chromosomes <- 5000
  chromosome_length <- 1
  recombination_rate <- 0
  mutation_rate <- 0
  smoothing_penalty <- seq(0, 30, 1)
}

dir.create(output_dir)

# ---------- (1) set up demographic model ---------- #
merger_time <- max_time * 1.0 #populations merge at this time (DISABLED)
migration_flip_time <- max_time * 0.4 #asymmetric migration rates flip at this time
migration_pulse_period <- 0 * max_time * c(0.15, 0.25) #migration from A->B pulses during this period (DISABLED)
migration_pulse <- 10^(-4.0) #migration rate from A->B during pulse
migration_high <- 10^(-5.0) #largest migration rate outside pulse
migration_low <- 10^(-6.0) #smallest of migration rate outside pulse
migration_growth_rate <- 0.001 #logistic rate of change for migration rate
average_ne <- 2e4 #Ne oscillates around this value
number_of_oscillations <- 3 #Ne shows this many oscillations before merger

time <- seq(0.0, max_time, length.out=1001)
start <- time[2:length(time)-1]
freq1 <- start*2*pi*number_of_oscillations/max_time
freq2 <- (start + 0.5*max_time/number_of_oscillations)*2*pi*number_of_oscillations/max_time
true_demographic_parameters <- array(0, c(2, 2, length(start)))
true_demographic_parameters[1,1,] <- ifelse(start > merger_time, average_ne, average_ne * (1 + 0.5*cos(freq1)))
true_demographic_parameters[2,2,] <- ifelse(start > merger_time, average_ne, average_ne * (1 + 0.5*cos(freq2)))
true_demographic_parameters[1,2,] <- ifelse(start > merger_time, 0, 
  migration_low + (migration_high - migration_low) * plogis(migration_growth_rate * (start - migration_flip_time)))
true_demographic_parameters[2,1,] <- ifelse(start > merger_time, 0, 
  migration_low + (migration_high - migration_low) * plogis(-migration_growth_rate * (start - migration_flip_time)))
true_demographic_parameters[1,2,] <- ifelse(start > migration_pulse_period[1] & start < migration_pulse_period[2],
                                            migration_pulse, true_demographic_parameters[1,2,])

population_tree <- PopulationTree$new(gsub("merger", merger_time, "(A:merger,B:merger);"), time_breaks=time)
population_tree$set_demographic_parameters(true_demographic_parameters)
save(population_tree, file=paste0(output_prefix, ".true_parameters.RData"))

# ---------- (2) simulate data ---------- #
set.seed(random_seed)

population_tree$msprime_simulate(
  outfile=output_prefix,
  sample_sizes=rep(haploids_per_population, 2), 
  chromosomes=number_of_chromosomes,
  chromosome_length=chromosome_length, 
  recombination_rate=recombination_rate,
  mutation_rate=mutation_rate,
  random_seed=sample.int(2^10, 1),
  what="tree_sequence"
)

population_tree$msprime_simulate(
  outfile=paste0(output_prefix, ".holdout"),
  sample_sizes=rep(haploids_per_population, 2), 
  chromosomes=holdout_chromosomes,
  chromosome_length=chromosome_length, 
  recombination_rate=recombination_rate,
  mutation_rate=mutation_rate,
  random_seed=sample.int(2^10, 1),
  what="tree_sequence"
)

# ---------- (3) calculate rates --------- #
# TODO: in the above simulated chromosomes are concatenated into a single .ts;
# but if reconstructed with relate or tsinfer will be separate .ts files; make
# example of merging ObservedTrioRates across chromosomes
epoch_breaks <- seq(0, max_time, length.out=number_of_epochs+1)
pop_assign <- rep(1:2-1, each=haploids_per_population)
sample_sets <- lapply(1:2-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- c("A", "B")

obs_rates <- ObservedTrioRates(
  ts = paste0(output_prefix, ".ts"),
  sample_sets = sample_sets, 
  time_breaks = epoch_breaks,
  bootstrap_blocks = 100
)
reticulate::py_save_object(obs_rates, filename=paste0(output_prefix, ".rates.pickle"))

obs_holdout_rates <- ObservedTrioRates(
  ts = paste0(output_prefix, ".holdout.ts"),
  sample_sets = sample_sets, 
  time_breaks = epoch_breaks,
  bootstrap_blocks = 100
)
reticulate::py_save_object(obs_holdout_rates, filename=paste0(output_prefix, ".holdout_rates.pickle"))

rates <- obs_rates$rates()
rates_sd <- obs_rates$std_dev(num_replicates=1000, random_seed=random_seed)
rates_denom <- obs_rates$denominator()[,1:ncol(rates)]
rates_precision <- smoothed_bootstrap_precision(
  mean=rates, 
  sd=rates_sd, 
  denominator=rates_denom, 
  method="gam"
)
rownames(rates) <- rownames(rates_precision) <- obs_rates$emission_states()
colnames(rates) <- colnames(rates_precision) <- obs_rates$epochs()

holdout_rates <- obs_holdout_rates$rates()
holdout_rates_sd <- obs_holdout_rates$std_dev(num_replicates=1000, random_seed=random_seed)
holdout_rates_denom <- obs_holdout_rates$denominator()[,1:ncol(holdout_rates)]
holdout_rates_precision <- smoothed_bootstrap_precision(
  mean=holdout_rates, 
  sd=holdout_rates_sd, 
  denominator=holdout_rates_denom, 
  method="gam"
)
rownames(holdout_rates) <- rownames(holdout_rates_precision) <- obs_holdout_rates$emission_states()
colnames(holdout_rates) <- colnames(holdout_rates_precision) <- obs_holdout_rates$epochs()

save(rates, rates_sd, rates_precision, rates_denom,
     holdout_rates, holdout_rates_sd, holdout_rates_precision, holdout_rates_denom,
     epoch_breaks, file=paste0(output_prefix, ".rates.RData"))

# ---------- (4) estimate demographic parameters --------- #
load(paste0(output_prefix, ".rates.RData"))

#if precision gets too large, gradient vanishes; so threshold.
#I should come up with a way to pick this automatically
rates_precision <- ifelse(rates_precision > 1e7, 1e7, rates_precision)
holdout_rates_precision <- ifelse(holdout_rates_precision > 1e7, 1e7, holdout_rates_precision)

upper_bound_migration <- 10^(-2) #if the migration rates/time unit get too high, matrix exponentiation will fail
lower_bound_migration <- 0
upper_bound_ne <- Inf
lower_bound_ne <- 10^2 #if the effective pop sizes/time unit get too low, matrix exponentiation will fail

population_tree <- PopulationTree$new(
  gsub("merger", merger_time, "(A:merger,B:merger);"), 
  time_breaks=epoch_breaks, 
  default_ne=1e4, #starting values for optimization
  default_migration_rate=5e-5 #starting values for optimization
)

lower_bounds <- array(log10(lower_bound_migration), c(2, 2, length(epoch_breaks)-1)) 
for(i in 1:2) lower_bounds[i,i,] <- log10(lower_bound_ne)
upper_bounds <- array(log10(upper_bound_migration), c(2, 2, length(epoch_breaks)-1)) 
for(i in 1:2) upper_bounds[i,i,] <- log10(upper_bound_ne)

fitted_models <- list()
for(pen in smoothing_penalty)
{
  cat("Penalty: ", pen, "\n", sep="")
  penalty <- matrix(pen, 2, 2) #smoothing penalty
  fitted_models[[as.character(pen)]] <- tryCatch({
    coaldecoder(
                coalescence_rates=rates,
                bootstrap_precision=rates_precision,
                epoch_durations=population_tree$epoch_durations(),
                demographic_parameters=population_tree$demographic_parameters(),
                admixture_coefficients=population_tree$admixture_coefficients(),
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                control=list(lmm=100),
                penalty=penalty,
                holdout=list(holdout_rates, holdout_rates_precision),
                debug_trace=TRUE,
                verbose=TRUE
                )
  }, error = function(e) {
    .coaldecoder_trace
  })
}

save(fitted_models, population_tree, file=paste0(output_prefix, ".fits.RData"))

## ----------- (6) figures ---------- #
load(paste0(output_prefix, ".true_parameters.RData"))
true_population_tree <- population_tree
load(paste0(output_prefix, ".fits.RData"))
load(paste0(output_prefix, ".rates.RData"))

library(ggplot2)
library(dplyr)
library(cowplot)
library(abind)

theme_custom <- function() 
{
  theme_cowplot() + 
    theme(text=element_text(family="Garamond"), 
          legend.text=element_text(size=8),
          legend.title=element_text(size=8),
          panel.border=element_rect(fill=NA),
          plot.title=element_text(size=14, face="plain"))
}

# (B) cross validation of smoothing penalty
cv_scores <- unlist(lapply(fitted_models, function(x) x$cross_validation_score))
best_penalty <- names(which.min(cv_scores))
cv_scores <- data.frame(score=cv_scores, penalty=as.numeric(names(cv_scores)))
cv_scores$weight <- exp(min(cv_scores$score) - cv_scores$score)
cv_scores %>% 
  ggplot() + 
  geom_line(aes(x=penalty, y=score)) +
  geom_point(aes(x=penalty, y=score, color=penalty, fill=penalty), pch=21) +
  theme_custom() +
  theme(legend.position="none") +
  scico::scale_fill_scico("Smoothing penalty", palette="vik", midpoint=as.numeric(best_penalty)) +
  scico::scale_color_scico("Smoothing penalty", palette="vik", midpoint=as.numeric(best_penalty)) +
  ggtitle("B. Cross-validation") +
  xlab("Smoothing penalty") + ylab("Cross-validation score") -> cv_score

# (A) fitted vs true parameter values, across range of smoothing penalty
parameters_to_df <- function(x, epoch_breaks)
{
  x <- abind(x, x[,,dim(x)[3]], along=3)
  coaldecoder:::.melt(log10(x)) %>%
    dplyr::rename(epoch=Var3, pop1=Var1, pop2=Var2) %>%
    dplyr::mutate(
      pop1 = c("A", "B")[pop1],
      pop2 = c("A", "B")[pop2],
      parameter = ifelse(pop1==pop2, 
        paste0("1/Size of ", pop1),
        paste0("Migr. ", pop1, "\U2192", pop2)
      ),
      type = ifelse(pop1==pop2, "Ne", "Migration"),
      time = epoch_breaks[epoch],
      value = ifelse(pop1==pop2, -value, value),
    ) %>%
    filter(!is.infinite(value))
}

parameters_df <- parameters_to_df(true_population_tree$demographic_parameters(), true_population_tree$.breaks)
parameters_df$penalty <- "truth"
for(i in names(fitted_models))
{
  tmp <- parameters_to_df(fitted_models[[i]]$demographic_parameters, epoch_breaks)
  tmp$penalty <- i
  parameters_df <- rbind(parameters_df, tmp)
}

time_scale = 1000
y_max <- 10^(-3.6)
x_max <- max(epoch_breaks)
parameters_df %>% filter(penalty=="truth") %>%
ggplot(aes(x=time/time_scale, group=parameter)) +
    geom_step(data=parameters_df %>% 
                filter(!(penalty %in% c('0', 'truth'))) %>% 
                mutate(penalty=as.numeric(penalty)),
              aes(y = 10^value, color=penalty, group=penalty), alpha=0.5) +
    geom_step(aes(y = 10^value), color="black", size=0.5) +
    scale_y_log10("Parameter value", breaks=10^seq(-10,10,1),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    xlim(0,x_max/time_scale) +
    xlab("Generations before present (thousands)") +
    ggtitle("A. True values vs estimates") +
    geom_text(data=parameters_df %>% filter(penalty=="truth" & epoch==1),
              aes(y = y_max, label = parameter),
              hjust=0, vjust=1, size=4, color="black", family="Garamond") + 
    scico::scale_color_scico("Smoothing penalty", palette="vik", 
                             limits=c(0,max(smoothing_penalty)),
                             breaks=c(0,as.numeric(best_penalty),max(smoothing_penalty)), 
                             midpoint=as.numeric(best_penalty)) +
    guides(color=guide_colorbar(title.position="top", 
                                frame.colour="black", 
                                ticks.colour="black",
                                barheight=unit(0.015,"npc"))) +
    facet_grid(pop1~pop2) +
    theme_custom() +
    theme(panel.grid=element_blank(),
          strip.text=element_blank(),
          legend.position=c(0.65,0.05),
          legend.direction="horizontal",
          text=element_text(family="Garamond"),
          strip.background=element_blank()) -> true_vs_fitted

# (C) fitted vs true vs observed rates
tmp0 <- parameters_df %>% filter(penalty==best_penalty)
tmp1 <- parameters_df %>% filter(penalty=='truth') %>% 
  filter(time %in% tmp0$time) %>% filter(time < max_time) %>%
  mutate(value = ifelse(pop1 == pop2, 10^(-value), 10^value))
#sample true parameters at beginning of epochs
true_parameters_coarse <- population_tree$demographic_parameters()
true_parameters_coarse[] <- tmp1$value
population_tree$set_demographic_parameters(true_parameters_coarse)
true_rates <- population_tree$expected_coalescence_rates()
save(population_tree, true_rates, file=paste0(output_prefix, ".true_rates.RData"))

rates_df <- data.frame(obs=c(fitted_models[[best_penalty]]$rates),
                       fitted=c(fitted_models[[best_penalty]]$fitted_rates),
                       true=c(true_rates))
ggplot(rates_df, aes(x=true)) +
  geom_point(aes(y=fitted), size=1, color="gray50") +
  scale_x_log10("True coalescence rate", breaks=10^seq(-10,10,1),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10("Fitted coalescence rate", breaks=10^seq(-10,10,1),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_abline(intercept=0, slope=1) +
  ggtitle("C. Fitted vs true rates") +
  theme_custom() -> fitted_rates


png(paste0(output_prefix, ".png"), height=6, width=10, units='in', res=300)
gridExtra::grid.arrange(true_vs_fitted, cv_score, 
                        fitted_rates, layout_matrix = matrix(c(1,1,2,3), nrow=2), widths=c(0.65,0.35))
dev.off()

