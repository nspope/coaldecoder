# ------------ test rate calculation with missing data ---------------------#
# If we mask exactly 1/2 of every tree, then rates should be identical to
# unmasked.
#
# python3 simulate_2pop_iid_missing.py --samples 10 --trees 100 --out missing_10x2.ts --seed 1024
# -------------------------------------------------------------------------#

library(coaldecoder)
reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

foo <- ObservedTrioRates("missing_10x2.ts", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,Inf))
bar <- ObservedTrioRates("missing_10x2.ts.mask", list("A"=0:2, "B"=3:5, "C"=6:8), c(0.0,Inf))

unmasked <- foo$y[,,1]/foo$n[,,1]
masked <- bar$y[,,1]/bar$n[,,1]

test_1.1 <- dplyr::near(unmasked, masked)

stopifnot(all(test_1.1[!is.na(test_1.1)]))

foo_boot <- foo$bootstrapped_rates(num_replicates=10, random_seed=1)
bar_boot <- bar$bootstrapped_rates(num_replicates=10, random_seed=1)

test_1.2 <- dplyr::near(foo_boot, bar_boot)

stopifnot(all(test_1.2))
