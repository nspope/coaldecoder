.melt <- function(df, value.name = "value") 
{
  df_long <- as.data.frame(as.table(df), responseName = value.name)
  cnames <- colnames(df_long)
  dnames <- dimnames(df)
  j <- 1
  for (i in cnames[grepl("^Var", cnames)]) {
    df_long[[i]] <- as.numeric(factor(df_long[[i]]))
    j <- j + 1
  }
  df_long
}

.plot_demographic_parameters <- function(
  parameters,
  epoch_durations = NULL,
  hessian = NULL,
  label_space = 0.2,
  time_scale = 1)
{
  library(ggplot2)
  library(dplyr)

  stopifnot(length(dim(parameters)) == 3)
  stopifnot(all(parameters >= 0))
  stopifnot(nrow(parameters) == ncol(parameters))

  if (is.null(rownames(parameters))) 
  {
    rownames(parameters) <- colnames(parameters) <- 1:nrow(parameters)-1
  }

  if (!is.null(epoch_durations))
  {
    stopifnot(length(epoch_durations) == dim(parameters)[3])
    stopifnot(all(epoch_durations > 0))
  } else {
    epoch_durations <- rep(1, dim(parameters)[3])
  }

  if (!is.null(hessian))
  {
    stopifnot(nrow(hessian) == length(parameters))
    std_err <- array(sqrt(diag(solve(hessian))), dim(parameters))
  } else { std_err <- array(NA, dim(parameters)) }

  epoch_start <- cumsum(c(0, epoch_durations))

  # to get last "step" to show up
  popnames <- rownames(parameters)
  parameters <- array(
    c(parameters, parameters[,,length(epoch_durations)]),
    dim=dim(parameters) + c(0,0,1)
  )
  rownames(parameters) <- colnames(parameters) <- popnames
  std_err <- array(
    c(std_err, std_err[,,length(epoch_durations)]),
    dim=dim(std_err) + c(0,0,1)
  )

  .melt(log10(parameters)) %>%
    dplyr::left_join(.melt(std_err, value.name = "stdErr"), 
                     by=c("Var1", "Var2", "Var3")) %>%
    dplyr::rename(epoch=Var3, pop1=Var1, pop2=Var2) %>% 
    dplyr::mutate(
      pop1_idx = pop1-1, 
      pop2_idx = pop2-1,
      pop1 = rownames(parameters)[pop1],
      pop2 = rownames(parameters)[pop2],
      parameter = ifelse(pop1_idx==pop2_idx, 
        paste0("1/N[", pop1_idx, "]"),
        paste0("M[", pop1_idx, ",", pop2_idx, "]")
      ),
      type = ifelse(pop1_idx==pop2_idx, "Ne", "Migration"),
      time = epoch_start[epoch],
      value = ifelse(pop1_idx==pop2_idx, -value, value),
    ) %>%
    group_by(parameter) %>%
    mutate(min_epoch = min(epoch[!is.infinite(value)])) %>%
    filter(epoch >= min_epoch) %>%
    filter(!is.infinite(value)) -> parameters_df

  parameters_df %>%
  ggplot(aes(x=time/time_scale, group=parameter)) +
    geom_vline(xintercept=0, lty='22', col="gray90") + 
    annotate(geom="text", y=Inf, x=0, label="Present ", vjust=0, hjust=1, angle=90, size=3, col="gray75") +
    geom_step(aes(y = 10^value)) +
    geom_point(data=parameters_df %>% filter(epoch == min_epoch),
               aes(y = 10^value, group=parameter)) +
    scale_y_log10("Parameter value", 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    xlim(-max(parameters_df$time/time_scale) * label_space, max(parameters_df$time/time_scale)) +
    xlab(if(time_scale == 1) "Time" else paste0("Time/", time_scale)) +
    geom_text(data=parameters_df %>% filter(epoch==min_epoch),
              aes(y = 10^value, label = parameter),
              hjust=1.2, vjust=0.5, size=3, color="black") + 
    facet_grid(pop1~pop2) +
    theme_bw() + 
    theme(panel.grid=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(face="bold")) -> demographic_plot

  if (!is.null(hessian))
  {
    demographic_plot <- demographic_plot +
      #geom_step(aes(y=10^(value-stdErr)), lty='22', size=0.5, alpha=0.5) +
      #geom_step(aes(y=10^(value+stdErr)), lty='22', size=0.5, alpha=0.5)
      pammtools::geom_stepribbon(aes(ymin=10^(value-stdErr), ymax=10^(value+stdErr)), alpha=0.15) 
  }

  demographic_plot
}

#TODO make rates optional and fitted required
.plot_expected_coalescence_rates <- function(
  rates,
  observed = NULL,
  epoch_durations = NULL,
  log_transform = FALSE,
  time_scale = 1)
{
  library(ggplot2)
  library(dplyr)

  if (!is.null(observed))
  {
    stopifnot(all(dim(rates) == dim(observed)))
  } else {
    observed <- matrix(NA, nrow(rates), ncol(rates))
  }
  rownames(observed) <- rownames(rates)
  colnames(observed) <- colnames(rates)

  if (!is.null(epoch_durations))
  {
    stopifnot(length(epoch_durations) == ncol(rates))
    stopifnot(all(epoch_durations > 0))
  } else {
    epoch_durations <- rep(1, ncol(rates))
  }

  epoch_start <- cumsum(c(0, epoch_durations))

  .melt(rates, value.name = "expected") %>%
    dplyr::left_join(.melt(observed, value.name = "observed")) %>% 
    dplyr::mutate(epoch=Var2, 
                  emission=rownames(rates)[Var1],
                  time=epoch_start[epoch],
                  observed=as.numeric(observed),
                  expected=as.numeric(expected)) %>%
    dplyr::mutate(observed=ifelse(is.na(observed) | observed == 0, NA, observed),
                  expected=ifelse(is.na(expected) | expected == 0, NA, expected)) %>%
    dplyr::mutate(Event=gsub("::.+", "", emission)=="t1",
                  Topology=gsub("t.::", "", emission)) %>%
    dplyr::mutate(Event=ifelse(Event, "First coalescence", "Second coalescence")) -> rates_df

    rates_df %>%
    ggplot(aes(x=time/time_scale)) +
    geom_line(aes(y=expected, color=Topology)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(face="bold")) +
    facet_grid(Event~., scales="free_y") +
    guides(color=guide_legend("Topology")) +
    xlab(ifelse(time_scale == 1, "Time", paste0("Time/", time_scale))) +
    ylab("Coalescence rate") -> output

  if (log_transform)
  {
    output <- output + 
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) 
  }

  if (!all(is.na(observed)))
  {
    output <- output +
      geom_point(aes(y=observed, color=Topology)) 
  }

  output
}
