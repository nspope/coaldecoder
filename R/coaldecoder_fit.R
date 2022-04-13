diagonal_precision <- function(x)
{
  #Matrix::Diagonal(ncol(x), x = 1/apply(x, 2, var))
  diag(1/apply(x, 2, var))
}

coaldecoder <- function(
  coalescence_rates,
  bootstrap_coalescence_rates,
  epoch_durations,
  demographic_parameters,
  admixture_coefficients = NULL,
  penalty = NULL,
  precision_matrix_estimator = diagonal_precision,
  lower_bounds = NULL,
  upper_bounds = NULL,
  max_restart = 1000,
  verbose = TRUE,
  calculate_hessian = FALSE)
{
  library(Matrix)

  if (is.null(admixture_coefficients))
  {
    admixture_coefficients <- array(0, dim(demographic_parameters))
    for (i in 1:nrow(admixture_coefficients))
    {
      admixture_coefficients[i,i,] <- 1.0
    }
  }

  if (is.null(penalty))
  {
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
    upper_bounds <- array(Inf, dim(demographic_parameters))
  }

  # check inputs
  decoder <- CoalescentDecoder$new(nrow(demographic_parameters), epoch_durations, TRUE)
  emission_names <- decoder$emission_states(rownames(demographic_parameters))

  # check bootstrap dims
  stopifnot(all(rownames(bootstrap_coalescence_rates) == rownames(coalescence_rates)))
  stopifnot(all(emission_names %in% rownames(coalescence_rates)))
  rate_mapping <- match(emission_names, rownames(coalescence_rates))
  coalescence_rates <- coalescence_rates[rate_mapping,,drop=FALSE]
  bootstrap_coalescence_rates <- bootstrap_coalescence_rates[rate_mapping,,,drop=FALSE]

  # build precision matrices
  precision_matrices <- lapply(1:length(epoch_durations),
    function(i) {
      # TODO: this is hacky. Better implement real checks.
      out <- matrix(0, length(emission_names), length(emission_names))
      boots <- t(bootstrap_coalescence_rates[,i,])
      valid <- !is.na(rowMeans(boots))
      stopifnot(sum(valid) > 3)
      boots <- boots[valid,]
      nonzero <- abs(colMeans(boots)) > nrow(boots)*.Machine$double.eps
      if (any(nonzero))
      {
        out[nonzero, nonzero] <- precision_matrix_estimator(boots[,nonzero])
      }
      invalid <- is.na(rowMeans(out)) | is.infinite(rowMeans(out))
      out[invalid,invalid] <- 0
      as(out, "dgCMatrix")
  })

  # optimization functions
  parameter_mapping <- which(.invert_diagonal(demographic_parameters) != 0)
  obj <- function(x) {
    demographic_parameters[parameter_mapping] <- 10^x
    lik <- decoder$loglikelihood(coalescence_rates, precision_matrices, 
             decoder$initial_state_vectors(), demographic_parameters, admixture_coefficients)
    prior <- decoder$smoothness_penalty(demographic_parameters, penalty, order=1)
    -(lik$loglikelihood + prior$penalty)
    }

  gra <- function(x) {
    demographic_parameters[parameter_mapping] <- 10^x
    lik <- decoder$loglikelihood(coalescence_rates, precision_matrices, 
             decoder$initial_state_vectors(), demographic_parameters, admixture_coefficients)
    prior <- decoder$smoothness_penalty(demographic_parameters, penalty, order=1)
    gradient <- -(lik$gradient$M + prior$gradient) * (demographic_parameters * log(10))
    gradient[parameter_mapping]
  }

  # L-BFGS-B with restarts
  fit <- optim(log10(demographic_parameters[parameter_mapping]), fn=obj, gr=gra, 
               lower=lower_bounds[parameter_mapping], upper=upper_bounds[parameter_mapping], method="L-BFGS-B")
  .skip_to_end <- FALSE
  for(i in 1:max_restart)
  {
    if (fit$convergence == 0) 
    { 
      break
    } else if (fit$message == "NEW_X") {
      fit <- tryCatch({
          optim(fit$par, fn=obj, gr=gra, method="L-BFGS-B",
                lower=lower_bounds[parameter_mapping], upper=upper_bounds[parameter_mapping])
      }, error = function(cond) {
          warning("Optimization failed")
          .skip_to_end <<- TRUE
          return(fit)
      })
      if (.skip_to_end) break
      if (verbose) {
        cat("L-BFGS-B restart:", i, "Deviance:", 2*fit$value, "\n")
      }
    } else {
      warning("Optimization failed")
      return(fit)
    }
  }

  if (verbose) cat(fit$message, "\n")

  demographic_parameters[parameter_mapping] <- 10^fit$par
  std_error <- array(NA, dim(demographic_parameters))

  # return hessian if desired
  if (calculate_hessian)
  {
    hessian <- numDeriv::jacobian(gra, fit$par)
    hessian <- (hessian + t(hessian))/2
    std_error[parameter_mapping] <- sqrt(diag(solve(hessian)))
  }

  list(demographic_parameters=demographic_parameters,
       std_error=std_error,
       loglikelihood=-fit$value,
       optimizer=fit,
       hessian=if(calculate_hessian) hessian else NULL)
}

#-----------------DEPRECATED API--------------------#

.coaldecoder_depr <- function(
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
  calculate_hessian = FALSE)
{
  library(Matrix)

  stopifnot(nrow(precision) == prod(dim(rates)))
  stopifnot(dim(demographic_parameters)[3] == length(epoch_duration))
  stopifnot(all(demographic_parameters >= 0))
  stopifnot(all(epoch_duration >= 0))
  stopifnot(all(rates >= 0))

  precision <- as(precision, "dgCMatrix")

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

  deco <- decoder_gaussian$new(nrow(demographic_parameters), 3, TRUE)
  obj <- function(mm) {
    mm <- array(10^mm, dim=dim(demographic_parameters))
    lik <- deco$deviance(deco$initial_states(), precision, rates, 
                         population_mapping, mm, epoch_duration)
    prior <- deco$random_walk_prior(mm, penalty, order=1)
    lik$deviance + prior$deviance
    }
  gra <- function(mm) {
    mm <- array(10^mm, dim=dim(demographic_parameters))
    lik <- deco$deviance(deco$initial_states(), precision, rates, 
                         population_mapping, mm, epoch_duration)
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
