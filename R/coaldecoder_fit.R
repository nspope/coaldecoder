coaldecoder <- function(
  rates,
  precision,
  demographic_parameters,
  epoch_duration,
  population_mapping = NULL,
  penalty = NULL,
  lower_bounds = NULL,
  upper_bounds = NULL,
  max_restart = 1000,
  verbose = TRUE,
  calculate_hessian = TRUE)
{
  stopifnot(nrow(precision) == prod(dim(rates)))
  stopifnot(dim(demographic_parameters)[3] == length(epoch_duration))
  stopifnot(all(demographic_parameters >= 0))
  stopifnot(all(epoch_duration >= 0))
  stopifnot(all(rates >= 0))

  if (class(precision) == "matrix")
  {
    precision <- as("dgCMatrix", precision)
  }

  if (!is.null(population_mapping))
  {
    stopifnot(nrow(demographic_parameters) == nrow(population_mapping))
    stopifnot(length(epoch_duration) == ncol(population_mapping))
  } else {
    population_mapping <- matrix(1:nrow(demographic_parameters)-1,
                                 nrow(demographic_parameters),
                                 length(epoch_duration))
  }

  if (!is.null(penalty))
  {
    stopifnot(nrow(penalty) == nrow(demographic_parameters))
    stopifnot(ncol(penalty) == ncol(demographic_parameters))
    stopifnot(all(penalty >= 0))
  } else {
    penalty <- matrix(0, nrow(demographic_parameters), ncol(demographic_parameters))
  }

  if (!is.null(lower_bounds))
  {
    stopifnot(all(dim(demographic_parameters) == dim(lower_bounds)))
  } else {
    lower_bounds <- array(-Inf, dim(demographic_parameters))
  }
  if (!is.null(upper_bounds))
  {
    stopifnot(all(dim(demographic_parameters) == dim(upper_bounds)))
  } else {
    upper_bounds <- array(-Inf, dim(demographic_parameters))
  }

  obj <- function(mm) {
    mm <- array(10^mm, dim=dim(demographic_parameters))
    lik <- deco$deviance(deco$initial_states(), precision, rates, remap, mm, epoch_duration)
    prior <- deco$random_walk_prior(mm, penalty, order=1)
    lik$deviance + prior$deviance
    }
  gra <- function(mm) {
    mm <- array(10^mm, dim=dim(demographic_parameters))
    lik <- deco$deviance(deco$initial_states(), precision, rates, remap, mm, epoch_duration)
    prior <- deco$random_walk_prior(mm, penalty, order=1)
    (lik$gradient + prior$gradient) * (mm * log(10))
  }
  fit <- optim(log10(demographic_parameters), fn=obj, gr=gra, method="L-BFGS-B")
  for(i in 1:max_restart)
  {
    if (fit$convergence == 0) 
    { 
      break
    } else if (fit$message == "NEW_X") {
      fit <- optim(fit$par, fn=obj, gr=gra, method="L-BFGS-B")
      if (verbose) {
        cat("L-BFGS-B restart:", i, "Deviance:", fit$value, "\n")
      }
    } else {
      warning("Optimization failed")
      return(fit)
    }
  }
  if (calculate_hessian)
  {
    hessian <- numDeriv::jacobian(gra, fit$par)
    fit$hessian <- (hessian + t(hessian))/2
  }
  fit
}
