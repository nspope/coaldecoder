
plot_demographic_model <- function(
  demographic_parameters,
  epoch_durations,
  hessian = NULL,
  type = c("both", "diagonal", "offdiagonal"),
  label_space = 0.2)
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

  melt <- function(df, value.name = "value") {
    as.data.frame(as.table(df), responseName = value.name) %>%
      dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Var")), list(factor)) %>%
      dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Var")), list(as.numeric)) %>% 
      as.data.frame
  }

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

  melt(log10(demographic_parameters)) %>%
    dplyr::left_join(melt(std_err, value.name = "stdErr"), 
                     by=c("Var1", "Var2", "Var3")) %>%
    dplyr::rename(epoch=Var3, pop1=Var1, pop2=Var2) %>% 
    dplyr::mutate(
      parameter = ifelse(pop1==pop2, 
        paste0("N[", pop1, "]"),
        paste0("M[", pop1, ",", pop2, "]")
      ),
      time = epoch_start[epoch]
    ) -> parameters_df

  if (type == "diagonal") {
    parameters_df <- parameters_df %>% filter(pop1 == pop2)
  } else if (type == "offdiagonal") {
    parameters_df <- parameters_df %>% filter(pop1 != pop2)
  }

  ggplot(parameters_df, aes(x=time, group=parameter)) +
    geom_step(aes(y = 10^value)) +
    scale_y_log10("Parameter value",
                  labels = scales::trans_format("log10", 
                                                scales::math_format(10^.x))
                  ) +
    xlim(0, max(parameters_df$time) * (1 + label_space)) +
    xlab("Time") +
    geom_text(data=parameters_df %>% filter(epoch==max(parameters_df$epoch)),
              aes(y = 10^value, 
                  label = parameter),
              hjust=0, size=3) + 
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
