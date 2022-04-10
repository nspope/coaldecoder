library(coaldecoder)

reticulate::source_python(system.file("python", "calculate_rates.py", package = "coaldecoder"))

#-----------no bootstrapping-------------#
#compare new rate calculator to old edge-diff version

{
foo <- ObservedTrioRates("example_3x3.ts", list(0:2, 3:5, 6:8), c(0.0,Inf))
old_rates <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5, 6:8), c(0.0,Inf))
new_rates <- foo$y/foo$n

test_1.1 <- all(dplyr::near(sort(unique(new_rates[1:18,1])), sort(unique(old_rates$y[old_rates$y > 0]))))
}

{
foo <- ObservedTrioRates("example_3x3.ts", list(0:2, 3:5, 6:8), c(0.0,10000,20000,30000))
old_rates <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5, 6:8), c(0.0,10000,20000,30000,Inf))
old_rates_nz <- calculate_rates(c("example_3x3.ts"), list(0:2, 3:5, 6:8), c(0.0,Inf))$y[,,1] > 0
old_rates <- apply(old_rates$y, 3, function(x) x[old_rates_nz])
old_rates <- apply(old_rates, 2, sort)
new_rates <- foo$y/foo$n
new_rates <- new_rates[1:(nrow(new_rates)/2),]
new_rates <- apply(new_rates, 2, sort)

test_1.2 <- all(dplyr::near(unlist(new_rates), unlist(old_rates)))
}

#-----------test bootstrapping-----------#
{
foo <- ObservedTrioRates("example_3x3.ts", list(0:2, 3:5, 6:8), c(0.0,10000,20000,30000), bootstrap_replicates=1, bootstrap_blocks=10, random_seed = 1)

test_2.1 <- all(dplyr::near(foo$n[,1], foo$n_boot[,1,1]))
test_2.2 <- !all(dplyr::near(foo$y[,1], foo$y_boot[,1,1]))
}

#-----------test merging-----------------#
{
foo <- ObservedTrioRates("example_3x3.ts", list(0:2, 3:5, 6:8), c(0.0,10000,20000,30000), bootstrap_replicates=1, bootstrap_blocks=10, random_seed = 1)
foo2 <- ObservedTrioRates("example_3x3.ts", list(0:2, 3:5, 6:8), c(0.0,10000,20000,30000), bootstrap_replicates=1, bootstrap_blocks=10, random_seed = 1)
foo2$y <- foo2$y * 0.75
foo2$n <- foo2$n * 0.75
foo2$y_boot <- foo2$y_boot * 0.75
foo2$n_boot <- foo2$n_boot * 0.75
foo$join(foo2)

test_3.1 <- all(dplyr::near( (foo2$y + foo2$y * 1/0.75)/2, foo$y))
test_3.2 <- all(dplyr::near( (foo2$n + foo2$n * 1/0.75)/2, foo$n))
test_3.3 <- all(dplyr::near( (foo2$y_boot + foo2$y_boot * 1/0.75)/2, foo$y_boot))
test_3.4 <- all(dplyr::near( (foo2$n_boot + foo2$n_boot * 1/0.75)/2, foo$n_boot))
}

#-----------test masking-----------------#
#TODO:

#-----------test ancient sample----------#
#TODO:
