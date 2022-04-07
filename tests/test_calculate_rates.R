library(coaldecoder)

#-----------no bootstrapping-------------#
#check equivalence with old method
old_rates <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5), c(0.0,1000,2000,3000))
rates <- calculate_trio_rates(c("example_3x3.ts"), list(0:2, 3:5), c(0.0,1000,2000,3000), merge_ts=TRUE)

deco <- CoalescentDecoder$new(2, diff(c(0.0,1000,2000,3000)), TRUE, TRUE)
rates <- deco$coalescence_rates(rates$y[,,1], rates$n[,1,drop=FALSE], diff(c(0.0,1000,2000,3000)), FALSE, FALSE)

or <- unique(c(old_rates$y))
nr <- unique(c(rates[1:6,]))

test_1.1 <- all(dplyr::near(or, nr))

#----------with bootstrapping, no merging--------------#
rates <- calculate_trio_rates(c("example_3x3.ts"), list(0:2, 3:5), c(0.0,1000,2000,3000), bootstrap_replicates=20, bootstrap_blocks=20, random_seed=1024, merge_ts=FALSE)
test_2.1 <- all(lapply(rates, function(x) class(x) == "TrioCoalescenceRates"))

#----------with bootstrapping, merging, input to decoder--------------#
rates <- calculate_trio_rates(c("example_3x3.ts"), list(0:2, 3:5), c(0.0,10000), bootstrap_replicates=10, bootstrap_blocks=10, random_seed=1024, merge_ts=TRUE)

test_3.1 <- class(rates) == "TrioCoalescenceRates"

deco <- CoalescentDecoder$new(2, rates$y, rates$n, diff(c(0.0,10000)), TRUE, TRUE)
