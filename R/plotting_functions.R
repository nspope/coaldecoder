.melt <- function(df, value.name = "value") 
{
  df_long <- as.data.frame(as.table(df), responseName = value.name)
  cnames <- colnames(df_long)
  dnames <- dimnames(df)
  j <- 1
  for (i in cnames[grepl("^Var", cnames)]) {
    df_long[[i]] <- as.numeric(factor(df_long[[i]]))
    if (!is.null(dnames[[j]]))
    {
      df_long[[i]] <- dnames[[j]][df_long[[i]]]
    }
    j <- j + 1
  }
  df_long
}

plot_demographic_model <- function(
  demographic_parameters,
  epoch_durations,
  hessian = NULL,
  type = c("both", "diagonal", "offdiagonal"),
  label_space = 0.2,
  time_scale = 1)
{
  type <- match.arg(type)

  library(ggplot2)
  library(dplyr)

  stopifnot(length(dim(demographic_parameters)) == 3)
  stopifnot(all(demographic_parameters >= 0))
  stopifnot(nrow(demographic_parameters) == ncol(demographic_parameters))

  if (!is.null(hessian))
  {
    stopifnot(nrow(hessian) == length(demographic_parameters))
    std_err <- array(sqrt(diag(solve(hessian))), dim(demographic_parameters))
  } else { std_err <- array(NA, dim(demographic_parameters)) }

  epoch_start <- cumsum(c(0, epoch_durations))

  # to get last "step" to show up
  demographic_parameters <- array(
    c(demographic_parameters, demographic_parameters[,,length(epoch_durations)]),
    dim=dim(demographic_parameters) + c(0,0,1)
  )
  std_err <- array(
    c(std_err, std_err[,,length(epoch_durations)]),
    dim=dim(std_err) + c(0,0,1)
  )

  .melt(log10(demographic_parameters)) %>%
    dplyr::left_join(.melt(std_err, value.name = "stdErr"), 
                     by=c("Var1", "Var2", "Var3")) %>%
    dplyr::rename(epoch=Var3, pop1=Var1, pop2=Var2) %>% 
    dplyr::mutate(
      pop1 = pop1-1, pop2=pop2-1,
      parameter = ifelse(pop1==pop2, 
        paste0("N[", pop1, "]"),
        paste0("M[", pop1, ",", pop2, "]")
      ),
      type = ifelse(pop1==pop2, "Ne", "Migration"),
      time = epoch_start[epoch]
    ) %>%
    group_by(parameter) %>%
    mutate(min_epoch = min(epoch[!is.infinite(value)])) %>%
    filter(epoch >= min_epoch) -> parameters_df

  if (type == "diagonal") {
    parameters_df <- parameters_df %>% filter(pop1 == pop2)
  } else if (type == "offdiagonal") {
    parameters_df <- parameters_df %>% filter(pop1 != pop2)
  }

  parameters_df %>%
  ggplot(aes(x=time/time_scale, group=parameter)) +
    geom_step(aes(y = 10^value)) +
    geom_point(data=parameters_df %>% filter(epoch == min_epoch),
               aes(y = 10^value, group=parameter)) +
    scale_y_log10("Parameter value", #sec.axis = sec_axis(~ (. - diag_scale$mean[1])/diag
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    #xlim(0, max(parameters_df$time) * (1 + label_space)) +
    xlim(-max(parameters_df$time/time_scale) * label_space, max(parameters_df$time/time_scale)) +
    xlab(if(time_scale == 1) "Time" else paste0("Time/", time_scale)) +
    geom_vline(xintercept=0, lty='22', col="gray75") + 
    annotate(geom="text", y=Inf, x=0, label="Present ", vjust=-0.5, hjust=1, angle=90, size=3, col="gray75") +
    geom_text(data=parameters_df %>% filter(epoch==min_epoch),
              aes(y = 10^value, label = parameter),
              hjust=1.2, vjust=0.5, size=3, color="black") + 
    facet_grid(pop1~pop2) +
    #facet_grid(type~., scales="free_y") +
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

plot_coalescence_rates <- function(
  rates,
  epoch_durations,
  fitted = NULL,
  log_transform = FALSE,
  time_scale = 1)
{
  library(ggplot2)
  library(dplyr)

  if (!is.null(fitted))
  {
    stopifnot(all(dim(rates) == dim(fitted)))
  } else {
    fitted <- matrix(NA, nrow(rates), ncol(rates))
  }
  rownames(fitted) <- rownames(rates)
  colnames(fitted) <- colnames(rates)

  epoch_start <- cumsum(c(0, epoch_durations))

  .melt(rates, value.name = "observed") %>%
    dplyr::left_join(.melt(fitted, value.name = "fitted")) %>%
    dplyr::mutate(epoch=as.numeric(Var2), 
                  time=epoch_start[epoch],
                  fitted=as.numeric(fitted),
                  observed=as.numeric(observed)) %>%
    dplyr::mutate(fitted=ifelse(is.na(fitted) | fitted == 0, NA, fitted),
                  observed=ifelse(is.na(observed) | observed == 0, NA, observed)) %>%
    dplyr::mutate(Topology=gsub("t.::", "", Var1),
                  Event=gsub("::.+", "", Var1)=="t1") %>%
    dplyr::mutate(Event=ifelse(Event, "First coalescence", "Second coalescence")) -> rates_df

    rates_df %>%
    ggplot(aes(x=time/time_scale)) +
    geom_point(aes(y=observed, color=Topology)) +
    geom_line(aes(y=fitted, color=Topology)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(face="bold")) +
    facet_grid(Event~., scales="free_y") +
    guides(color=guide_legend("Topology")) +
    xlab(if(time_scale == 1) "Time" else paste0("Time/", time_scale)) +
    ylab("Coalescence rate") -> plot

  if (log_transform)
  {
    plot <- plot + 
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) 
  }

  plot
}
