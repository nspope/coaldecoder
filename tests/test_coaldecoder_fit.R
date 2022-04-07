library(coaldecoder)

#------------test 2-population fit with no penalty--------#
#we're just checking to see if this runs without failing

rates <- calculate_trio_rates(c("example_3x3.ts"), list(0:2, 3:5), c(0.0,10000,20000), bootstrap_replicates=20, bootstrap_blocks=20, random_seed=1024, merge_ts=TRUE)

M_start <- array(NA, c(2,2,2))
M_start[1,1,] <- 10000
M_start[2,2,] <- 10000
M_start[1,2,] <- 1e-5
M_start[2,1,] <- 1e-5

A <- array(0, c(2,2,2))
A[1,1,] <- 1
A[2,2,] <- 1

foo <- coaldecoder_trio(rates, M_start, A, diff(c(0.0,10000,20000)))



