.invert_diagonal <- function(arr)
{
  stopifnot(length(dim(arr)) == 3)
  for (i in 1:dim(arr)[3])
  {
    diag(arr[,,i]) <- 1/diag(arr[,,i])
  }
  arr
}

PopulationTree <- setRefClass("PopulationTree", fields = c(".prefix", ".M", ".A", ".tree", ".breaks", ".mask", ".times"))

PopulationTree$methods(
initialize = function(topology, time_breaks, default_ne = 1e4, default_migration_rate = 1e-4) 
{ 
  .prefix <<- "[PopulationTree] "

  ##important to test: ancient populations that merge immediately!
  ##simultaneous merger
  #time_breaks <- seq(0,10000,1000)
  #topology <- c("((A:2000,B:950):2000,C:4000):3000;")
  #default_migration <- 1e-4
  #default_ne <- 1e4

  #what if I try to fit a model where the sampling time is incorrectly specified between observed/demographic?
  #should I save .population_times in both and check these? probably

  stopifnot(is.character(topology))
  if (substr(topology, nchar(topology), nchar(topology)) != ";")
  {
    topology <- paste0(topology, ";")
  }

  stopifnot(is.numeric(time_breaks))
  stopifnot(!is.unsorted(time_breaks))
  stopifnot(all(time_breaks >= 0) & (0 %in% time_breaks))

  tree <- ape::reorder.phylo(ape::read.tree(text = topology), "postorder")

  stopifnot(!duplicated(tree$tip.label))

  node_depth <- ape::node.depth.edgelength(tree)
  tree$root.edge <- max(max(time_breaks) - max(node_depth), 0)
  node_depth <- max(node_depth) - node_depth
  stopifnot(all(!is.na(node_depth)))
  population_times <- node_depth[1:length(tree$tip.label)]
  names(population_times) <- tree$tip.label

  # snap internal nodes to time_breaks
  node_dist <- outer(node_depth, time_breaks[-length(time_breaks)], "-")
  node_epoch <- rep(NA, length(node_depth))
  node_is_internal <- rep(c(FALSE,TRUE), c(length(tree$tip.label), tree$Nnode))
  for (i in 1:nrow(node_dist))
  {
    times <- node_dist[i,]
    if (0 %in% times) { 
      node_epoch[i] <- which(times == 0)
    } else if (all(times > 0)) {
      node_epoch[i] <- NA
      warning(paste0(.prefix, "Node occurs after oldest epoch, ignoring"))
    } else {
      if (node_is_internal[i])
      {
        node_epoch[i] <- which.min(abs(times))
        warning(paste0(.prefix, "Node does not coincide with any epoch boundary, snapping to nearest epoch"))
      } else {
        stop(paste0(.prefix, "Population sampling times must coincide with an epoch boundary"))
      }
    }
  }

  # identify epochs where there are mergers, and who mergers into who
  node_labels <- c(tree$tip.label, rep(NA, tree$Nnode))
  merger <- list()
  for (parent in unique(tree$edge[,1]))
  {
    children <- tree$edge[tree$edge[,1]==parent, 2]
    node_labels[parent] <- sort(node_labels[children])[1]
    merger[[parent]] <- sort(node_labels[children])
  }
  visit_order <- order(node_epoch)

  # make parameter arrays
  # internally we store Ne as 1/Ne
  M <- array(default_migration_rate, c(length(tree$tip.label), length(tree$tip.label), length(time_breaks) - 1))
  A <- array(0, c(length(tree$tip.label), length(tree$tip.label), length(time_breaks) - 1))
  rownames(M) <- colnames(M) <- rownames(A) <- colnames(A) <- tree$tip.label
  dimnames(M)[[3]] <- dimnames(A)[[3]] <- 
    paste0("[", time_breaks[2:length(time_breaks)-1], ",", time_breaks[2:length(time_breaks)], ")")
  for(i in 1:nrow(M)) { M[i,i,] <- 1/default_ne; A[i,i,] <- 1 }

  for (i in visit_order)
  {
    if (i <= length(tree$tip.label)) #is a tip
    {
      if (node_epoch[i] > 1) #is sampled before the present day
      {
        idx <- match(node_labels[i], tree$tip.label)
        M[idx,-idx,2:node_epoch[i]-1] <- 0
        M[-idx,idx,2:node_epoch[i]-1] <- 0
        M[idx,idx,2:node_epoch[i]-1] <- 0
      }
    } else {
      if (!is.na(node_epoch[i])) #merger occurs within span of time_breaks
      {
        idx <- match(merger[[i]][-c(1)], tree$tip.label)
        M[idx,-idx,node_epoch[i]:dim(M)[3]] <- 0
        M[-idx,idx,node_epoch[i]:dim(M)[3]] <- 0
        M[idx,idx,node_epoch[i]:dim(M)[3]] <- 0

        jdx <- match(merger[[i]][1], tree$tip.label)
        A[jdx,idx,node_epoch[i]] <- 1
        A[idx,idx,node_epoch[i]] <- 0
      }
    }
  }

  M <- M[order(rownames(M)), order(colnames(M)), , drop=FALSE]
  A <- A[order(rownames(A)), order(colnames(A)), , drop=FALSE]
  population_times <- population_times[order(names(population_times))]

  .M <<- M
  .A <<- A
  .tree <<- tree
  .breaks <<- time_breaks
  .mask <<- M == 0
  .times <<- population_times
})

#------------------ Getters ---------------#
PopulationTree$methods(demographic_parameters = function() { return(.invert_diagonal(.M)) })
PopulationTree$methods(admixture_coefficients = function() { return(.A) })
PopulationTree$methods(time_breaks = function() { return(.breaks) })
PopulationTree$methods(epoch_durations = function() { return(diff(.breaks)) })
PopulationTree$methods(population_sampling_times = function() { return(.times) })
PopulationTree$methods(population_tree = function() { return(.tree) })
PopulationTree$methods(population_names = function() { return(names(.times)) })

#------------------ Setters ---------------#
PopulationTree$methods(set_demographic_parameters = function(demographic_parameters)
{
  stopifnot(length(dim(demographic_parameters)) == 3)
  stopifnot(all(dim(demographic_parameters) == dim(.M)))
  stopifnot(all(!is.na(demographic_parameters)))

  demographic_parameters <- .invert_diagonal(demographic_parameters)

  stopifnot(all(is.finite(demographic_parameters))) #fails if input population size is 0
  if (any((demographic_parameters == 0) != .mask)) warning("Masking post-merger demographic parameters")

  demographic_parameters[.mask] <- 0
  dimnames(demographic_parameters) <- dimnames(.M)

  .M <<- demographic_parameters
})

#----------------- Interface to coaldecoder ------------------#

PopulationTree$methods(expected_coalescence_rates = function()
{
  decoder <- CoalescentDecoder$new(length(population_names()), epoch_durations(), TRUE)
  rates <- decoder$expected_rates(decoder$initial_state_vectors(), demographic_parameters(), admixture_coefficients())$y
  rownames(rates) <- decoder$emission_states(population_names())
  colnames(rates) <- paste0("(", time_breaks()[2:length(time_breaks())-1], 
                            ",", time_breaks()[2:length(time_breaks())], "]")
  rates
})

PopulationTree$methods(occupancy_probabilities = function()
{
  decoder <- CoalescentDecoder$new(length(population_names()), epoch_durations(), TRUE)
  occupancy <- decoder$occupancy_probabilities(demographic_parameters(), admixture_coefficients())
  rownames(occupancy) <- rownames(demographic_parameters())
  colnames(occupancy) <- colnames(demographic_parameters())
  occupancy
})

#------------------ Interface to msprime ---------------------#
PopulationTree$methods(msprime_simulate = function(outfile, sample_sizes, chromosomes=1, chromosome_length=1e6, recombination_rate=1e-8, mutation_rate=0.0, what=c("tree_sequence", "msprime_inputs"), random_seed=NULL)
{
  what <- match.arg(what)

  reticulate::source_python(system.file("python", "simulate.py", package = "coaldecoder"))

  msprime_simulate (sample_sizes, population_sampling_times(), population_names(),
                    demographic_parameters(), admixture_coefficients(), 
                    epoch_durations(), chromosome_length, recombination_rate,
                    chromosomes, outfile, what, mutation_rate, random_seed)
})


#------------------ Plotting functions -----------------------#
PopulationTree$methods(plot_population_tree = function()
{
  plot(.tree, root.edge=TRUE)
  tree_depth <- max(ape::node.depth.edgelength(.tree))
  offset <- max(tree_depth - max(.breaks), 0)
  axis(1, at=offset+.breaks, line=1, labels=FALSE)
  axis(1, at=range(offset+.breaks), line=1, tick=FALSE, labels=range(.breaks)[2:1])
})

PopulationTree$methods(plot_expected_coalescence_rates = function(...)
{
  .plot_expected_coalescence_rates(rates=expected_coalescence_rates(), epoch_durations=epoch_durations(), ...)
})

PopulationTree$methods(plot_demographic_parameters = function(...)
{
  .plot_demographic_parameters(parameters=demographic_parameters(), epoch_durations=epoch_durations(), ...)
})

PopulationTree$methods(plot_occupancy_probabilities = function(...)
{
  .plot_occupancy_probabilities(occupancy=occupancy_probabilities(), epoch_durations=epoch_durations(), ...)
})
