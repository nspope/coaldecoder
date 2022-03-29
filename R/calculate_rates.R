
calculate_rates <- function(ts_paths,
                            sample_sets,
                            time_breaks,
                            bootstrap_replicates = 0,
                            bootstrap_blocks = 1,
                            random_seed = NULL,
                            merge_ts = FALSE)
{
  parse_raw_inputs <- function(raw)
  {
    deco <- decoder_gaussian$new(length(sample_sets), 3, TRUE)
    states <- deco$emission_states()
    classes <- deco$emission_classes()

    y <- array(0, c(nrow(states), nrow(classes), length(time_breaks)-1))
    n <- array(0, c(nrow(classes), 1))
    y_boot <- array(0, c(dim(y), dim(raw$y_boot)[3]))
    n_boot <- array(0, c(dim(n), 1, dim(raw$y_boot)[3]))

    for(i in 1:nrow(raw$labels))
    {
      x <- raw$labels[i,]

      ec <- sapply(0:max(raw$labels), function(n) sum(x==n))
      ec <- which(apply(classes, 1, function(y) all(y==ec)))
      es <- sort(x[1:2])
      es <- which(apply(states, 1, function(y) all(y==es)))
      stopifnot(length(es)==1)
      stopifnot(length(ec)==1)

      y[es,ec,] <- y[es,ec,] + raw$y[i,]
      n[ec] <- n[ec] + raw$n[i]

      if (dim(raw$y_boot)[3] > 0) 
      {
        for (j in 1:dim(raw$y_boot)[3])
        {
          y_boot[es,ec,,j] <- y_boot[es,ec,,j] + raw$y_boot[i,,j]
          n_boot[ec,1,j] <- n_boot[ec,1,j] + raw$n_boot[i,j]
        }
      }
    }
    list(y=y, 
         n=n, 
         y_boot=y_boot, 
         n_boot=n_boot, 
         file=raw$path, 
         trees=raw$trees, 
         size=raw$size, 
         classes=paste0("{", apply(classes, 1, paste, collapse=","), "}"),
         states=paste0(apply(classes, 1, paste, collapse="x")),
         NULL
    )
  } 

  reticulate::source_python(system.file("python", "coaldecoder.py", package = "coaldecoder"))
  rates <- lapply(
    trio_first_coalescence_rates(ts_paths, sample_sets, time_breaks, 
                                 bootstrap_replicates, bootstrap_blocks, 
                                 random_seed),
    parse_raw_inputs)

  if (merge_ts)
  {
    rates <- Reduce(function(a,b)
    {
      stopifnot(all(a$classes == b$classes))
      stopifnot(all(a$states == b$states))
      size <- a$size + b$size
      a$y <- a$y * a$size/size + b$y * b$size/size
      a$n <- a$n * a$size/size + b$n * b$size/size
      a$y_boot <- a$y_boot * a$size/size + b$y_boot * b$size/size
      a$n_boot <- a$n_boot * a$size/size + b$n_boot * b$size/size
      a$trees <- a$trees + b$trees
      a$file <- c(a$file, b$file)
      a$size <- size
    }, rates)
  }

  rates
}
