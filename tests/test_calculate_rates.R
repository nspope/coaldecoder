library(coaldecoder)

reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

#-----------no bootstrapping-------------#
#compare new rate calculator to old edge-diff version

test_against_old <- function()
{
foo <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,Inf))
old_rates <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5, 6:8), c(0.0,Inf))
new_rates <- foo$y[,,1]/foo$n[,,1]

test_1.1 <- all(dplyr::near(sort(unique(new_rates[1:18,1])), sort(unique(old_rates$y[old_rates$y > 0]))))

stopifnot(test_1.1)
}

test_against_old_windowed <- function()
{
foo <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000))
old_rates <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5, 6:8), c(0.0,10000,20000,30000,Inf))
old_rates_nz <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5, 6:8), c(0.0,Inf))$y[,,1] > 0
old_rates <- apply(old_rates$y, 3, function(x) x[old_rates_nz])
old_rates <- apply(old_rates, 2, sort)
new_rates <- foo$y[,,1]/foo$n[,,1]
new_rates <- new_rates[1:(nrow(new_rates)/2),]
new_rates <- apply(new_rates, 2, sort)

test_1.2 <- all(dplyr::near(unlist(new_rates), unlist(old_rates)))

stopifnot(test_1.2)
}

#-----------test bootstrapping-----------#
test_bootstrapping <- function()
{
foo <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000), bootstrap_blocks = 10)
boot_y <- foo$bootstrapped_rates(num_replicates=1000, random_seed=1)
y <- foo$rates()
y_hat <- apply(boot_y, c(1,2), mean)
plot(y, y_hat); abline(0,1) #good

test_2.1 <- all( abs(y - y_hat) < 1e-5, na.rm=TRUE)

std_d <- foo$std_dev(num_replicates=1000, random_seed=1)
boot_r <- foo$bootstrapped_rates(num_replicates=1000, random_seed=1)
boot_sd <- apply(boot_r, c(1,2), sd)

test_2.2 <- all(dplyr::near(std_d, boot_sd), na.rm=TRUE)

stopifnot(test_2.1)
stopifnot(test_2.2)
}

#-----------test block-by-trees----------#
test_blocking <- function()
{
foo <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000), trees_per_block = 50)

stopifnot(foo$bootstrap_blocks == floor(foo$num_trees / 50))
}

#-----------test merging-----------------#
test_merging <- function()
{
foo <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000))
foo2 <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000), bootstrap_blocks=10)
foo3 <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000), bootstrap_blocks=5)

foo2$join(foo3)

test_3.1 <- sum(foo2$block_span) == 1
test_3.2 <- length(foo2$block_span) == 15
test_3.3 <- foo2$bootstrap_blocks == length(foo2$block_span)
test_3.4 <- foo2$sequence_length == 2*foo$sequence_length
test_3.5 <- all(dplyr::near(foo$rates(), foo2$rates()), na.rm=TRUE)

stopifnot(test_3.1)
stopifnot(test_3.2)
stopifnot(test_3.3)
stopifnot(test_3.4)
stopifnot(test_3.5)
}

#-----------test precision smoothing-----#
test_smoothing <- function()
{
#this isnt a test, just playing around with ideas for implementation
# python3 simulate_3pop_iid_ancient.py --samples 10 --trees 50000 --out asympt_iid_ancient_10x3.ts --seed 1024
filename <- "asympt_iid_ancient_10x3.ts"
num_pops <- 3
haps_per_pop <- 10
time_breaks_true <- c(as.matrix(read.table(paste0(filename, ".t"))))
time_breaks <- seq(0.0, max(time_breaks_true), 500)
              
pop_assign <- rep(1:num_pops-1, each=haps_per_pop)
sample_sets <- lapply(1:num_pops-1, function(i) as.integer(which(pop_assign == i)-1)) 
names(sample_sets) <- LETTERS[1:num_pops]
obs <- ObservedTrioRates(filename, sample_sets, time_breaks, trees_per_block=100)
ee <- obs$rates()
oo <- obs$std_dev(1000, 1)
aa <- apply(obs$n[,-ncol(obs$n),], c(1,2), sum)

df <- melt(log10(oo))
df$n <- c(aa)
df$x <- c(log10(ee))
df <- df[is.finite(rowSums(df)),]
df$Var1 <- factor(df$Var1)
library(gamm4)
library(mgcv)
oof0 <- mgcv::gam(value ~ x + s(Var1,n,bs='fs',m=1), data=df) #works well...
oof1 <- mgcv::gamm(value ~ x + s(Var1,n,bs='fs',m=1), data=df) 
oof2 <- gamm4::gamm4(value ~ x + s(Var1,n,bs='fs',m=1), data=df)
df$fit0 <- predict(oof0)
df$fit1 <- predict(oof1$gam)
df$fit2 <- predict(oof2$mer)
plot(df$value, df$fit0)
plot(df$value, df$fit1)
plot(df$value, df$fit2)

library(ggplot2)
library(dplyr)
library(reshape2)
dff <- melt(log10(oo))
dff$dur <- rep(diff(time_breaks),each=nrow(oo))
dff$x <- c(log10(ee))
dff$n <- c(aa)
dff %>% group_by(Var1) %>% mutate(nn = n/max(n)) %>%
  ggplot() + geom_line(aes(x=nn, y=value,color=factor(Var1)))

#thoughts: use a GAM to smooth std deviation, with rate, denominator as covariates?

}

#-----------test masking-----------------#
test_masking <- function()
{
ts_obj <- reticulate::py_run_string("
import tskit
import numpy as np
ts = tskit.load('example_3x3.ts')
breaks = np.array([i for i in ts.breakpoints()])
", local=TRUE)
breaks <- ts_obj[["breaks"]]

mask <- cbind(breaks[2:length(breaks) - 1], breaks[2:length(breaks) - 1] + 1)
foo <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000), mask=mask)

test_1 <- dplyr::near(foo$sequence_length, ts_obj[["ts"]][["sequence_length"]] - nrow(mask))
stopifnot(test_1)

mask <- cbind(breaks[2:length(breaks) - 1], breaks[2:length(breaks)])
mask <- mask[1:10,]
bar <- ObservedTrioRates("example_3x3.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,10000,20000,30000), mask=mask)

test_2 <- dplyr::near(bar$sequence_length, max(breaks) - breaks[11])
test_3 <- bar$num_trees == length(breaks) - 1 - 10
stopifnot(test_2)
stopifnot(test_3)
}

#-----------test ancient sample----------#
#TODO:

#----------------------------------------#
test_against_old()
test_against_old_windowed()
test_bootstrapping()
test_blocking()
test_merging()
test_masking()

