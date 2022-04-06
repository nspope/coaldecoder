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
