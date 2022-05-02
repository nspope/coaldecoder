smoothed_bootstrap_precision <- function(mean, sd)
{
  stopifnot(length(dim(mean)) == 2)
  stopifnot(length(dim(sd)) == 2)
  stopifnot(all(dim(mean) == dim(sd)))

  df <- data.frame(x=log(c(mean)), y=log(c(sd)))
  valid <- is.finite(df$x) & is.finite(df$y)
  fit <- lm(y~x, data=df[valid,])
  precision <- sd
  precision[valid] <- exp(-fitted(fit))
  precision
}

coaldecoder <- function(
  coalescence_rates,
  bootstrap_precision,
  epoch_durations,
  demographic_parameters,
  admixture_coefficients = NULL,
  penalty = NULL,
  lower_bounds = NULL,
  upper_bounds = NULL,
  state = NULL,
  max_restart = 100,
  verbose = TRUE,
  calculate_hessian = FALSE,
  holdout = NULL,
  debug_trace = FALSE)
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
  #TODO: check input types for rates, bootstrap precision
  stopifnot(all(dim(bootstrap_precision) == dim(coalescence_rates)))
  stopifnot(all(rownames(bootstrap_precision) == rownames(coalescence_rates)))
  stopifnot(all(emission_names %in% rownames(coalescence_rates)))
  rate_mapping <- match(emission_names, rownames(coalescence_rates))
  coalescence_rates <- coalescence_rates[rate_mapping,,drop=FALSE]
  bootstrap_precision <- bootstrap_precision[rate_mapping,,drop=FALSE]

  # calculate sqrt(precision) from bootstraps
  #bootstrap_precision <- 1./apply(bootstrap_coalescence_rates, c(1,2), sd)
  bootstrap_precision <- ifelse(is.na(bootstrap_precision) | is.infinite(bootstrap_precision), 0, bootstrap_precision)

  # optimization functions
  parameter_mapping <- which(.invert_diagonal(demographic_parameters) > 0.0)

  if (is.null(state))
  {
    state <- decoder$initial_state_vectors()
  } else {
    stopifnot(all(dim(state) == dim(decoder$initial_state_vectors())))
  }

  obj <- function(x) {
    demographic_parameters[parameter_mapping] <- 10^x
    if (debug_trace) #deleteme eventually
    {
      .coaldecoder_trace <<- list(demographic_parameters=demographic_parameters, 
                                  admixture_coefficients=admixture_coefficients,
                                  rates=coalescence_rates, 
                                  precision=bootstrap_precision, 
                                  state=state,
                                  epoch_durations=epoch_durations)
    }
    lik <- decoder$loglikelihood(coalescence_rates, bootstrap_precision,
             state, demographic_parameters, admixture_coefficients)
    prior <- decoder$smoothness_penalty(demographic_parameters, penalty, order=1)
    gradient <- -(lik$gradient$M + prior$gradient) * (demographic_parameters * log(10))
    .coaldecoder_gradient <<- gradient[parameter_mapping]
    return (-lik$loglikelihood - prior$penalty)
  }

  gra <- function(x) {
    #demographic_parameters[parameter_mapping] <- 10^x
    #lik <- decoder$loglikelihood(coalescence_rates, bootstrap_precision,
    #         decoder$initial_state_vectors(), demographic_parameters, admixture_coefficients)
    #prior <- decoder$smoothness_penalty(demographic_parameters, penalty, order=1)
    #gradient <- -(lik$gradient$M + prior$gradient) * (demographic_parameters * log(10))
    #gradient[parameter_mapping]
    return (.coaldecoder_gradient)
  }

  #TODO: use a more modern implementation of L-BFGS-B
  # L-BFGS-B with restarts
  fit <- optim(log10(demographic_parameters[parameter_mapping]), fn=obj, gr=gra, 
               lower=lower_bounds[parameter_mapping], upper=upper_bounds[parameter_mapping], method="L-BFGS-B")
  .coaldecoder_terminate <<- FALSE
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
          .coaldecoder_terminate <<- TRUE
          return(fit)
      })
      if (.coaldecoder_terminate) break
      if (verbose) {
        cat("L-BFGS-B restart:", i, "Deviance:", 2*fit$value, "\n")
      }
    } else {
      warning("Optimization failed")
      break
    }
  }

  if (verbose) cat(fit$message, "\n")

  demographic_parameters[parameter_mapping] <- 10^fit$par
  std_error <- array(NA, dim(demographic_parameters))

  # return hessian if desired
  if (calculate_hessian)
  {
    #TODO: test that scoping works here, could default to commented-out part of 'gra' above
    hessian <- numDeriv::jacobian(function(x) {obj(x); gra(x)}, fit$par)
    hessian <- (hessian + t(hessian))/2
    try({
      std_error[parameter_mapping] <- sqrt(diag(solve(hessian)))
    })
  }

  # cross-validation score
  cross_validation_score <- NA
  if (!is.null(holdout))
  {
    #TODO: clean this up
    rates_holdout <- holdout[[1]]
    bootstrap_holdout <- holdout[[2]]
    stopifnot(all(rownames(bootstrap_holdout) == rownames(rates_holdout)))
    stopifnot(all(emission_names %in% rownames(rates_holdout)))
    rate_mapping <- match(emission_names, rownames(rates_holdout))
    rates_holdout <- rates_holdout[rate_mapping,,drop=FALSE]
    bootstrap_holdout <- bootstrap_holdout[rate_mapping,,,drop=FALSE]
    
    bootstrap_precision_holdout <- 1./apply(bootstrap_holdout, c(1,2), sd)
    bootstrap_precision_holdout <- ifelse(is.na(bootstrap_precision_holdout) | is.infinite(bootstrap_precision_holdout), 0, bootstrap_precision_holdout)
    lik_holdout <- decoder$loglikelihood(rates_holdout, bootstrap_precision_holdout,
      state, demographic_parameters, admixture_coefficients)
    cross_validation_score <- -2 * lik_holdout$loglikelihood
  }

  decoder_fit <- decoder$expected_rates(state, demographic_parameters, admixture_coefficients)
  fitted_rates <- decoder_fit$y
  rownames(fitted_rates) <- rownames(coalescence_rates)
  colnames(fitted_rates) <- colnames(coalescence_rates)

  list(demographic_parameters=demographic_parameters,
       std_error=std_error,
       loglikelihood=-fit$value,
       optimizer=fit,
       fitted_rates=fitted_rates,
       rates=coalescence_rates,
       precision=bootstrap_precision,
       state=decoder_fit$X,
       hessian=if(calculate_hessian) hessian else NULL,
       cross_validation_score=if(!is.null(holdout)) cross_validation_score else NULL)
}

#optimize one epoch at a time
.coaldecoder_piecewise <- function(
  coalescence_rates,
  bootstrap_precision,
  epoch_durations,
  demographic_parameters,
  admixture_coefficients = NULL,
  state = NULL,
  penalty = NULL,
  lower_bounds = NULL,
  upper_bounds = NULL,
  numerical_gradient = FALSE,
  max_restart = 100,
  verbose = TRUE)
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
  #TODO: check input types for rates, bootstrap precision
  stopifnot(all(dim(bootstrap_precision) == dim(coalescence_rates)))
  stopifnot(all(rownames(bootstrap_precision) == rownames(coalescence_rates)))
  stopifnot(all(emission_names %in% rownames(coalescence_rates)))
  rate_mapping <- match(emission_names, rownames(coalescence_rates))
  coalescence_rates <- coalescence_rates[rate_mapping,,drop=FALSE]
  bootstrap_precision <- bootstrap_precision[rate_mapping,,drop=FALSE]

  # calculate sqrt(precision) from bootstraps
  #bootstrap_precision <- 1./apply(bootstrap_coalescence_rates, c(1,2), sd)
  bootstrap_precision <- ifelse(is.na(bootstrap_precision) | is.infinite(bootstrap_precision), 0, bootstrap_precision)

  # optimization functions
  epoch_fits <- list()
  loglik <- list()

  fitted_rates <- matrix(NA, nrow(coalescence_rates), ncol(coalescence_rates))
  rownames(fitted_rates) <- rownames(coalescence_rates)

  if (is.null(state))
  {
    state <- decoder$initial_state_vectors()
  } else {
    stopifnot(all(dim(state) == dim(decoder$initial_state_vectors())))
  }

  last <- NULL

  for (i in 1:dim(demographic_parameters)[3])
  {
    if (verbose) cat("Epoch ", i, "\n", sep="")

    ldecoder <- CoalescentDecoder$new(nrow(demographic_parameters), epoch_durations[i], TRUE)
    pars <- demographic_parameters[,,i,drop=FALSE]
    admix <- admixture_coefficients[,,i,drop=FALSE]
    rates <- coalescence_rates[,i,drop=FALSE]
    prec <- bootstrap_precision[,i,drop=FALSE]
    lb <- lower_bounds[,,i,drop=FALSE]
    ub <- upper_bounds[,,i,drop=FALSE]

    parameter_mapping <- which(.invert_diagonal(pars) > 0.0)

    obj <- function(x) {
      pars[parameter_mapping] <- 10^x
      lik <- ldecoder$loglikelihood(rates, prec, state, pars, admix)
      prior <- list(penalty=0, gradient=array(0,dim(pars)))
      if (!is.null(last))
      {
        prior$penalty <- sum((-0.5 * c(penalty)^2 * c(log10(pars) - log10(last))^2)[parameter_mapping])
        prior$gradient <- array(c(penalty)^2 * c(log10(pars) - log10(last)), dim(pars))
      }
      gradient <- -(lik$gradient$M) * pars * log(10) - prior$gradient
      .coaldecoder_gradient <<- gradient[parameter_mapping]
      #..store <<- list(rates=rates, prec=prec, state=state, pars=pars, admix=admix, deco=ldecoder)
      return (-lik$loglikelihood - prior$penalty)
    }

    gra <- function(x) {
      return (.coaldecoder_gradient)
    }

    .coaldecoder_terminate <<- FALSE
    fit <- optim(
            log10(pars[parameter_mapping]), fn=obj,
            gr=if(numerical_gradient) NULL else gra,
            lower=lb[parameter_mapping], upper=ub[parameter_mapping],
            method="L-BFGS-B"
    )
    for(j in 1:max_restart)
    {
      if (fit$convergence == 0) 
      { 
        break
      } else if (fit$message == "NEW_X") {
        fit <- tryCatch({
          optim(fit$par, fn=obj, 
                gr=if(numerical_gradient) NULL else gra, 
                lower=lb[parameter_mapping], upper=ub[parameter_mapping],
                method="L-BFGS-B")
        }, error = function(cond) {
          warning("Optimization failed")
          .coaldecoder_terminate <<- TRUE
          return(fit)
        })
        if (.coaldecoder_terminate) break
        if (verbose) {
          cat("L-BFGS-B restart:", j, "Deviance:", 2*fit$value, "\n")
        }
      } else {
        warning("Optimization failed")
        break
      }
    }

    if (verbose) cat(fit$message, "\n")

    pars[parameter_mapping] <- 10^fit$par

    ##<DEBUG>
    #browser()
    #foo <- numDeriv::grad(function(x) obj(x), log10(pars))
    #bar <- obj(log10(pars))
    #bar <- gra()
    #print(foo)
    #print(bar)
    #print("--")
    ##</DEBUG>

    last <- pars
    demographic_parameters[,,i] <- pars
    loglik[[i]] <- ldecoder$loglikelihood(rates, prec, state, pars, admix)
    fitted_rates[,i] <- loglik[[i]]$y_hat
    state <- loglik[[i]]$X
    epoch_fits[[i]] <- fit
  }

  list(demographic_parameters=demographic_parameters,
       fitted_rates=fitted_rates,
       rates=coalescence_rates,
       precision=bootstrap_precision,
       state=state,
       loglikelihood=loglik,
       optimizer=epoch_fits)
}

.coaldecoder_piecewise_admixture <- function(
  coalescence_rates,
  bootstrap_precision,
  epoch_durations,
  demographic_parameters,
  admixture_coefficients = NULL,
  state = NULL,
  penalty = NULL,
  lower_bounds = NULL,
  upper_bounds = NULL,
  numerical_gradient = FALSE,
  max_restart = 100,
  verbose = TRUE)
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
  #TODO: check input types for rates, bootstrap precision
  stopifnot(all(dim(bootstrap_precision) == dim(coalescence_rates)))
  stopifnot(all(rownames(bootstrap_precision) == rownames(coalescence_rates)))
  stopifnot(all(emission_names %in% rownames(coalescence_rates)))
  rate_mapping <- match(emission_names, rownames(coalescence_rates))
  coalescence_rates <- coalescence_rates[rate_mapping,,drop=FALSE]
  bootstrap_precision <- bootstrap_precision[rate_mapping,,drop=FALSE]

  # calculate sqrt(precision) from bootstraps
  #bootstrap_precision <- 1./apply(bootstrap_coalescence_rates, c(1,2), sd)
  bootstrap_precision <- ifelse(is.na(bootstrap_precision) | is.infinite(bootstrap_precision), 0, bootstrap_precision)

  # optimization functions
  epoch_fits <- list()
  loglik <- list()

  fitted_rates <- matrix(NA, nrow(coalescence_rates), ncol(coalescence_rates))
  rownames(fitted_rates) <- rownames(coalescence_rates)

  if (is.null(state))
  {
    state <- decoder$initial_state_vectors()
  }

  for (i in 1:dim(demographic_parameters)[3])
  {
    if (verbose) cat("Epoch ", i, "\n", sep="")

    ldecoder <- CoalescentDecoder$new(nrow(demographic_parameters), epoch_durations[i], TRUE)
    pars <- demographic_parameters[,,i,drop=FALSE]
    admix <- admixture_coefficients[,,i,drop=FALSE]
    rates <- coalescence_rates[,i,drop=FALSE]
    prec <- bootstrap_precision[,i,drop=FALSE]
    lb <- lower_bounds[,,i,drop=FALSE]
    ub <- upper_bounds[,,i,drop=FALSE]

    parameter_mapping <- which(.invert_diagonal(pars) > 0.0)
    parameter_mapping2 <- which(lower.tri(admix[,,1]) | upper.tri(admix[,,1]))

    col_softmax <- function(z) {
      for (j in 1:ncol(z)) {
        z[,j] <- exp(z[,j]) / sum(exp(z[,j]))
      }
      z
    }

    jac_softmax <- function(z, y) {
      d <- matrix(0, nrow(y), ncol(y))
      for (j in 1:ncol(z))
      {
        jac <- -outer(z[,j], z[,j])
        diag(jac) <- z[,j] * (1-z[,j])
        d[,j] <- jac %*% y[,j]
      }
      d
    }

    #TODO: this will optimize all admixture parameters
    admix_template <- admix*0
    admix_start <- diag(nrow(admix))
    for(j in 1:nrow(admix)) {
      admix_start[j,j] <- 0.9
      admix_start[-j,j] <- 0.1 / nrow(admix)
    }
    lb2 <- matrix(-Inf, nrow(admix_start), ncol(admix_start))
    ub2 <- matrix(Inf, nrow(admix_start), ncol(admix_start))

    obj <- function(x) {
      p1 <- 1:length(parameter_mapping)
      p2 <- (length(parameter_mapping)+1):(length(parameter_mapping)+length(parameter_mapping2))
      pars[parameter_mapping] <- 10^x[p1]
      admix_template[parameter_mapping2] <- x[p2]
      admix_template[,,1] <- col_softmax(admix_template[,,1])
      lik <- ldecoder$loglikelihood(rates, prec, state, pars, admix_template)
      gradient <- -lik$gradient$M * pars * log(10)
      gradient2 <- -jac_softmax(admix_template[,,1], lik$gradient$A[,,1])
      .coaldecoder_gradient <<- c(gradient[parameter_mapping], gradient2[parameter_mapping2])
      return (-lik$loglikelihood)
    }

    gra <- function(x) {
      return (.coaldecoder_gradient)
    }

    .coaldecoder_terminate <<- FALSE
    fit <- optim(
            c(log10(pars[parameter_mapping]),log(admix_start[parameter_mapping2])),
            fn=obj,
            gr=if(numerical_gradient) NULL else gra,
            lower=c(lb[parameter_mapping], lb2[parameter_mapping2]), 
            upper=c(ub[parameter_mapping], ub2[parameter_mapping2]),
            method="L-BFGS-B"
    )
    for(j in 1:max_restart)
    {
      if (fit$convergence == 0) 
      { 
        break
      } else if (fit$message == "NEW_X") {
        fit <- tryCatch({
            optim(
            fit$par,
            fn=obj,
            gr=if(numerical_gradient) NULL else gra,
            lower=c(lb[parameter_mapping], lb2[parameter_mapping2]), 
            upper=c(ub[parameter_mapping], ub2[parameter_mapping2]),
            method="L-BFGS-B"
            )
        }, error = function(cond) {
          warning("Optimization failed")
          .coaldecoder_terminate <<- TRUE
          return(fit)
        })
        if (.coaldecoder_terminate) break
        if (verbose) {
          cat("L-BFGS-B restart:", j, "Deviance:", 2*fit$value, "\n")
        }
      } else {
        warning("Optimization failed")
        break
      }
    }

    if (verbose) cat(fit$message, "\n")

    p1 <- 1:length(parameter_mapping)
    p2 <- (length(parameter_mapping)+1):(length(parameter_mapping)+length(parameter_mapping2))
    pars[parameter_mapping] <- 10^fit$par[p1]
    demographic_parameters[,,i] <- pars
    admix_template[parameter_mapping2] <- fit$par[p2]
    admix[,,1] <- col_softmax(admix_template[,,1])
    admixture_coefficients[,,i] <- admix[,,1]
    loglik[[i]] <- ldecoder$loglikelihood(rates, prec, state, pars, admix)
    fitted_rates[,i] <- loglik[[i]]$y_hat
    state <- loglik[[i]]$X
    epoch_fits[[i]] <- fit
  }

  list(demographic_parameters=demographic_parameters,
       admixture_coefficients=admixture_coefficients,
       fitted_rates=fitted_rates,
       rates=coalescence_rates,
       precision=bootstrap_precision,
       state=state,
       loglikelihood=loglik,
       optimizer=epoch_fits)
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
