## requires profiling to be enabled:
## uncomment lines in Makevars.in, profiler.cpp and rebuild

  library(coaldecoder)
  library(Matrix)

  M <- matrix(0.5, 31, 31)

  start_profiler("profile_constructor.prof")
  hmm <- CoalescentDecoder$new(nrow(M), c(1000), TRUE)
  stop_profiler()

  start_profiler("profile_rate_matrix_0.prof")
  Q1 <- hmm$transition_rates(M)
  stop_profiler()

  start_profiler("profile_rate_matrix_1.prof")
  Q2 <- hmm$transition_rates_devel(M, FALSE)
  stop_profiler()

  start_profiler("profile_rate_matrix_2.prof")
  Q3 <- hmm$transition_rates_devel(M, TRUE)
  stop_profiler()

  all(dplyr::near(Q1@x, Q3@x))
  all(dplyr::near(Q1@x, Q2@x))

  ## diagonal admixture is very quick, not worth optimizing
  A <- diag(31)
  start_profiler("profile_admixture_matrix_2.prof")
  A0 <- hmm$admixture_proportions(A)
  stop_profiler()

  system.time(hmm$admixture_proportions(A))
