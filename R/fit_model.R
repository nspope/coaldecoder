map_vector_to_cube <- function(x, y, mapping)
{
  stopifnot(all(dim(mapping) == dim(y)))
  stopifnot(length(unique(mapping[!is.na(mapping)])) == length(x))
  y[!is.na(mapping)] <- x[mapping[!is.na(mapping)]]
  y
}

map_cube_to_vector <- function(y, f, mapping)
{
  stopifnot(all(dim(mapping) == dim(y)))
  ag <- aggregate(y[!is.na(mapping)], by=list(ind = mapping[!is.na(mapping)]), f)
  x <- rep(NA, nrow(ag))
  x[ag$ind] <- ag$x
  x
}

make_parameter_mapping <- function(model_matrix, asymmetric = TRUE, fix.ne = FALSE)
{
  pars <- c()
  mapping <- array(NA, dim(model_matrix))
  for(t in 1:dim(model_matrix)[3])
  {
    offset <- length(pars)
    active <- which(lower.tri(model_matrix[,,t], diag = !fix.ne) | (asymmetric & upper.tri(model_matrix[,,t])))
    new_pars <- offset + 1:length(active)
    mapping[,,t][active] <- new_pars
    if (!asymmetric){
      mapping[,,t][upper.tri(mapping[,,t])] <- t(mapping[,,t])[upper.tri(mapping[,,t])]
    } 
    pars <- c(pars, new_pars)
  }
  mapping
}

fit_model <- function(n, y, remap, model_matrix, duration, penalty = c(0, 0, 0, 0, 0), penalty.order = 2, ne.bounds = c(3, 7), m.bounds = c(-Inf, -2), asymmetric = TRUE, fix.ne = FALSE, mapping = NULL, max.eval = 5000, conv.tol = 1e-12, use.approx.rate.matrix = TRUE, verbose = 100, num_ind = 3)
{
  num_pop <- dim(model_matrix)[1]
  #num_ind <- floor(-1/2 + sqrt(1/4 + 2*nrow(y)))#wrong this gives # pop
  deco <- decoder$new(num_pop, num_ind, use.approx.rate.matrix)

  mapping_ <- if(is.null(mapping)) make_parameter_mapping(model_matrix, asymmetric = asymmetric, fix.ne = fix.ne) else mapping
  x0_ <- log10(model_matrix)
  state_ <- NA
  verbose_ <- verbose
  gradient_ <- NA
  penalty_ <- penalty
  order_ <- penalty.order
  lb_ = array(apply(model_matrix, 3, function(x) { x[] <- m.bounds[1]; diag(x) <- ne.bounds[1]; x }), dim(model_matrix))
  ub_ = array(apply(model_matrix, 3, function(x) { x[] <- m.bounds[2]; diag(x) <- ne.bounds[2]; x }), dim(model_matrix))
  iter_ <- 0
  opts_ <- list(algorithm="NLOPT_LD_LBFGS", maxeval=max.eval, xtol_rel=0., xtol_abs=0., ftol_abs=conv.tol, ftol_rel=conv.tol) 
  print(system.time(initial_ <- deco$deviance(deco$initial_states(), n, y, remap, 10^x0_, duration)))

  fit <- nloptr::nloptr(
    eval_f = function(x, initial_states, size, counts, remap, model_matrix, duration) 
    {
      model_matrix <- map_vector_to_cube(10^x, model_matrix, mapping_)
      loss <-
        deco$deviance(initial_states, size, counts, remap, model_matrix, duration)[c("deviance","gradient")]
      gradient <- loss[["gradient"]]
      deviance <- loss <- loss[["deviance"]]
      #--smoothness penalty
      migr_fun <- deco$migration_function_penalty(diag(nrow(model_matrix)), remap, model_matrix, duration, penalty_, order_)
      loss <- loss + migr_fun[["deviance"]]
      gradient <- gradient + migr_fun[["gradient"]]
      #--snip
      #state_ <<- model_matrix
      gradient_ <<- gradient * model_matrix * log(10)
      iter_ <<- iter_ + 1
      if (verbose_ && iter_ %% verbose_ == 0) cat(iter_, deviance, loss, "\n")
      loss
    },
    eval_grad_f = function(x, initial_states, size, counts, remap, model_matrix, duration) 
    {
      gradient <- map_cube_to_vector(gradient_, sum, mapping_)
      gradient_ <<- NA #ensure function gets called beforehand
      gradient
    },
    opts = opts_,
    x0 = map_cube_to_vector(x0_, mean, mapping_),
    lb = map_cube_to_vector(lb_, min, mapping_),
    ub = map_cube_to_vector(ub_, max, mapping_),
    initial_states = deco$initial_states(),
    size = n,
    counts = y,
    remap = remap,
    model_matrix = model_matrix,
    duration = duration
    )
  
  model_matrix <- map_vector_to_cube(10^fit$solution, model_matrix, mapping_)

  convergence <- fit$message

  fit <- deco$deviance(deco$initial_states(), n, y, remap, model_matrix, duration)
  fit$migration_function <- deco$migration_function_penalty(diag(nrow(model_matrix)), remap, model_matrix, duration, penalty_, order_)

  attr(fit, "duration") <- duration
  attr(fit, "x0") <- x0_
  attr(fit, "model_matrix") <- model_matrix
  attr(fit, "mapping") <- mapping_
  attr(fit, "remap") <- remap
  attr(fit, "penalty") <- penalty
  attr(fit, "penalty.order") <- penalty.order
  attr(fit, "y") <- y
  attr(fit, "n") <- n
  attr(fit, "convergence") <- convergence
  attr(fit, "num_ind") <- num_ind
  attr(fit, "use.approx.rate.matrix") <- use.approx.rate.matrix

  if (verbose) print(convergence)
  fit
}

force_population_merger_1 <- function(model_matrix, mapping, epoch, from, to, ub = 0.01, lb = 1e-12)
{
  #force all lineages to immediately transition to an absorbing population
  model_matrix[from,from,epoch:dim(model_matrix)[3]] <- model_matrix[to,to,epoch:dim(model_matrix)[3]]#force ne to absorbing population
  model_matrix[from,to,epoch:dim(model_matrix)[3]] <- ub #force all lineages in "from" to move to "to" backwards in time
  model_matrix[from,-c(to,from),epoch:dim(model_matrix)[3]] <- lb #no migration of lineages to other demes aside from "to"
  model_matrix[-c(from),from,epoch:dim(model_matrix)[3]] <- lb #no migration of lineages from other demes into "from"

  new_mapping <- mapping

  new_mapping[from,from,epoch:dim(model_matrix)[3]] <- NA #ne in merged deme is nonidentifiable 
  new_mapping[from,to,epoch:dim(model_matrix)[3]] <- NA #force all lineages in "from" to move to "to" backwards in time
  new_mapping[from,-c(to,from),epoch:dim(model_matrix)[3]] <- NA #no migration of lineages to other demes aside from "to"
  new_mapping[-c(from),from,epoch:dim(model_matrix)[3]] <- NA #no migration of lineages from other demes into "from"

  new_par_index <- unique(new_mapping[!is.na(new_mapping)])
  new_mapping[!is.na(new_mapping)] <- match(new_mapping[!is.na(new_mapping)], new_par_index)

  list(model_matrix = model_matrix, mapping = new_mapping)
}

force_population_merger_2 <- function(model_matrix, mapping, epoch, from, to, lb = 1e-12)
{
  #isolate population in entrance direction
  model_matrix[-c(from),from,epoch:dim(model_matrix)[3]] <- lb #no migration of lineages from other demes into "from"

  new_mapping <- mapping

  new_mapping[-c(from),from,epoch:dim(model_matrix)[3]] <- NA #no migration of lineages from other demes into "from"

  new_par_index <- unique(new_mapping[!is.na(new_mapping)])
  new_mapping[!is.na(new_mapping)] <- match(new_mapping[!is.na(new_mapping)], new_par_index)

  list(model_matrix = model_matrix, mapping = new_mapping)
}

force_population_merger_3 <- function(model_matrix, mapping, epoch, from, to, lb = 1e-12)
{
  #isolate population in both directions
  model_matrix[from,-c(from),epoch:dim(model_matrix)[3]] <- lb #exit rates
  model_matrix[-c(from),from,epoch:dim(model_matrix)[3]] <- lb #entrance rates

  new_mapping <- mapping
  new_mapping[from,-c(from),epoch:dim(model_matrix)[3]] <- NA 
  new_mapping[-c(from),from,epoch:dim(model_matrix)[3]] <- NA 

  new_par_index <- unique(new_mapping[!is.na(new_mapping)])
  new_mapping[!is.na(new_mapping)] <- match(new_mapping[!is.na(new_mapping)], new_par_index)

  list(model_matrix = model_matrix, mapping = new_mapping)
}

plot_system <- function(fit, true_model = NULL, what = 2)
{
 library(ggplot2)
 library(reshape2)
 library(dplyr)

 fit_model_matrix <- attr(fit, "model_matrix")
 if (is.null(true_model))
 {
   model_matrix <- 10^attr(fit, "x0") 
 } else {
   model_matrix <- true_model
 }
 num_pop <- nrow(model_matrix)
 num_ind <- attr(fit, "num_ind")
 duration <- attr(fit, "duration")
 use.approx.rate.matrix <- attr(fit, "use.approx.rate.matrix")

 deco <- decoder$new(num_pop, num_ind, use.approx.rate.matrix)

 if (what == 1) {

 melt(log10(model_matrix), value.name = "true") %>% 
   dplyr::left_join(melt(log10(fit_model_matrix), value.name = "fit"), by=c("Var1","Var2","Var3")) %>%
   ggplot(aes(x=Var3, group=factor(Var2), color=factor(Var2))) + 
   geom_line(aes(y=fit)) +
   geom_line(aes(y=true),linetype=2) +
   facet_wrap(~Var1, ncol=num_pop, nrow=num_pop, scales="free_y") +
   theme_bw() 

 } else if (what == 2) {

 melt(deco$migration_function(model_matrix, duration), value.name = "true") %>% 
   dplyr::left_join(melt(deco$migration_function(fit_model_matrix, duration), value.name = "fit"), by=c("Var1","Var2","Var3")) %>%
   ggplot(aes(x=Var3,group=factor(Var1),color=factor(Var1))) + 
   geom_line(aes(y=fit)) +
   geom_line(aes(y=true),linetype=2) +
   facet_wrap(~Var2, ncol=num_pop, nrow=num_pop)

 } else if (what == 3) {

 melt(log10(model_matrix), value.name = "true") %>% 
   dplyr::left_join(melt(log10(fit_model_matrix), value.name = "fit"), by=c("Var1","Var2","Var3")) %>%
   dplyr::filter(Var1 == Var2) %>% 
   ggplot(aes(x=Var3, group=factor(Var2), color=factor(Var2))) + 
   geom_line(aes(y=fit)) +
   geom_line(aes(y=true),linetype=2)

 } else if (what == 4) {

 melt(deco$migration_function(fit_model_matrix, duration), value.name = "fit") %>%
   dplyr::group_by(Var1, Var2) %>%
   dplyr::arrange(Var3) %>%
   dplyr::summarise(Var3 = 2:length(Var3), dfit = diff(fit)) %>% 
   ggplot(aes(x=Var3,group=factor(Var1),color=factor(Var1))) + 
   geom_line(aes(y=dfit)) +
   facet_wrap(~Var2, ncol=num_pop, nrow=num_pop) 

 } else if (what == 5) {

 melt(deco$migration_operator(fit_model_matrix, duration), value.name = "fit") %>%
   dplyr::filter(Var1 != Var2) %>%
   ggplot(aes(x=Var3,group=factor(Var1),color=factor(Var1))) + 
   geom_line(aes(y=fit)) +
   facet_wrap(~Var2, ncol=num_pop, nrow=num_pop)

 } else if (what == 6) {

 melt(array(apply(deco$migration_operator(fit_model_matrix, duration),3,solve),dim(fit_model_matrix)), value.name = "fit") %>%
   dplyr::filter(Var1 != Var2) %>%
   ggplot(aes(x=Var3,group=factor(Var1),color=factor(Var1))) + 
   geom_line(aes(y=fit)) +
   facet_wrap(~Var2, ncol=num_pop, nrow=num_pop)
 }
}

plot_data <- function(model_fit, pairs, migr_window=NA, migr_spike=NA)
{

parse_state_names <- function(st, cl, nm)
{
  colnames(cl) <- nm
  out <- data.frame()
  for(i in 1:nrow(st))
    for(j in 1:nrow(cl))
    {
      a=nm[st[i,]+1]
      b=cl[j,]
      for(k in a) b[k] = b[k] - 1
      if(any(b<0)) {x <- NA} else if(sum(b>1)) {x <- NA} else {x <- paste0("((", paste(a,collapse=","), "),", names(which(b>0)), ")")}
      y=paste0("{", paste(sort(rep(colnames(cl), cl[j,])), collapse=","), "}")
      out <- rbind(out, data.frame(estate=i, eclass=j, tree=x, config=y))
    }
  out
}

duration <- attr(model_fit, "duration")
n <- attr(model_fit, "n")
y <- attr(model_fit, "y")
remap <- attr(model_fit, "remap")
model_matrix <- attr(model_fit, "model_matrix")
num_ind <- attr(model_fit, "num_ind")

if(!is.na(migr_window))
{
  model_matrix[pairs[1],pairs[2],migr_window] <- migr_spike
}

deco <- decoder$new(nrow(model_matrix),num_ind,TRUE)

fit <- deco$deviance(deco$initial_states(), n, y, remap, model_matrix, duration)

emission <- array(apply(fit$emission, 3, function(x) prop.table(x,2)),dim(fit$emission))
emission <- emission[-nrow(emission),,]
predicted <- fit$predicted
predicted <- predicted[-nrow(predicted),,]

#labels
emission_labels <- parse_state_names(deco$emission_states(), deco$emission_classes(), 1:nrow(model_matrix))

library(dplyr)
library(reshape2)
library(ggplot2)

melt(emission, value.name = "observed") %>% 
  left_join(melt(predicted, value.name = "fit")) %>%
  mutate(estate=Var1, eclass=Var2, time=Var3) %>% 
  select(estate, eclass, time, observed, fit) %>%
  left_join(emission_labels) %>%
  filter(!is.na(tree)) -> obs_vs_fit

keep <- grepl(pairs[1], obs_vs_fit$tree) | grepl(pairs[2], obs_vs_fit$tree)
not_pairs <- (1:nrow(model_matrix))[!((1:nrow(model_matrix)) %in% pairs)]
drop <- rep(FALSE, nrow(obs_vs_fit))
for(i in not_pairs) drop <- (drop | grepl(i, obs_vs_fit$tree))

obs_vs_fit %>% filter(!drop & keep) %>%
  ggplot(aes(x=time, color=tree)) +
    geom_point(aes(y=observed), size=1, alpha=1) +
    geom_line(aes(y=fit), linetype=1, position=position_dodge(width=0.2)) +
    theme_bw() +
    facet_wrap(~tree, scales="free") +
    theme(panel.grid=element_blank(), text=element_text(family="Garamond")) -> one

melt(log10(model_matrix), value.name = "parameters") -> parameters

rbind(parameters) %>%
  filter(Var1 %in% pairs & Var2 %in% pairs) %>%
  mutate(group = paste0(Var1,Var2)) %>%
  ggplot(aes(x = Var3, group=group, linetype=group)) + 
    geom_step(aes(y = parameters)) +
    theme_bw() +
    facet_wrap(~group, scales="free") +
    theme(panel.grid=element_blank(), text=element_text(family="Garamond")) -> two

  list(one, two)
}

add_ghost <- function(n, y, num_ind, num_pop)
{
  deco1 <- decoder$new(num_pop+1, num_ind, TRUE)
  deco2 <- decoder$new(num_pop, num_ind, TRUE)
  states_1 <- deco1$emission_states()
  states_2 <- deco2$emission_states()
  rows <- match(apply(states_2,1,paste,collapse="_"), apply(states_1,1,paste,collapse="_"))
  class_1 <- deco1$emission_classes()
  class_1 <- class_1[,1:num_pop]
  class_2 <- deco2$emission_classes()
  cols <- match(apply(class_2,1,paste,collapse="_"), apply(class_1,1,paste,collapse="_"))
  stopifnot(all(!is.na(rows)))
  stopifnot(all(!is.na(cols)))
  #
  x0 <- array(0, c(num_pop+1,num_pop+1,dim(y)[3]))
  x0 <- array(apply(x0, 3, function(x) { x[] <- -4; diag(x) <- 5; x }), dim(x0))
  remap <- matrix(1:(num_pop+1)-1, num_pop+1, dim(y)[3])
  #
  y_new=array(0,c(nrow(states_1),nrow(class_1),dim(y)[3])) 
  n_new=matrix(0,nrow(class_1),1)
  y_new[rows,cols,] <- y
  n_new[cols] <- n
  list(y=y_new,n=n_new,x0=x0,remap=remap)
}

subset_data <- function(n, y, num_ind, num_pop, new_pop, zero_which = NA)
{
  new_pop <- sort(new_pop)
  deco1 <- decoder$new(num_pop, num_ind, TRUE)
  deco2 <- decoder$new(length(new_pop), num_ind, TRUE)
  states_1 <- deco1$emission_states()
  states_2 <- deco2$emission_states()
  states_2 <- cbind(new_pop[states_2[,1]+1],new_pop[states_2[,2]+1])
  rows <- match(apply(states_2,1,paste,collapse="_"), apply(states_1,1,paste,collapse="_"))
  class_1 <- deco1$emission_classes()
  class_1 <- class_1[,new_pop+1]
  class_2 <- deco2$emission_classes()
  cols <- match(apply(class_2,1,paste,collapse="_"), apply(class_1,1,paste,collapse="_"))
  stopifnot(all(!is.na(rows)))
  stopifnot(all(!is.na(cols)))
  x0 <- array(0, c(length(new_pop),length(new_pop),dim(y)[3]))
  x0 <- array(apply(x0, 3, function(x) { x[] <- -4; diag(x) <- 5; x }), dim(x0))
  remap <- matrix(1:length(new_pop)-1,length(new_pop),dim(y)[3])
  #
  y=y[rows,cols,,drop=FALSE]
  n=n[cols,drop=FALSE]
  if (!is.na(zero_which))
  {
    zero_class <- class_2[,zero_which] > 0
    y[,zero_class,] <- 0
    n[zero_class] <- 0
  }
  list(y=y,n=n,labels=new_pop,x0=x0,remap=remap)
}

subset_and_refit <- function(fit, new_pop, zero_which = NA, ...)
{
  foo <- subset_data(attr(fit,"n"), attr(fit,"y"), attr(fit,"num_ind"), nrow(attr(fit,"model_matrix")), new_pop, zero_which)
  fit_model(foo$n, foo$y, foo$remap, 10^foo$x0, attr(fit,"duration"), attr(fit,"penalty"), ...)
}
