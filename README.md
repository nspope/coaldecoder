# coaldecoder

Demographic inference from reconstructed genealogies using trio first coalescence rates. Note that this package is at a 'proof of concept' stage, the API and implementation are bound to change.

# Python dependencies:

- tskit (>0.4.0)
- numpy

# R dependencies:

- RcppArmadillo
- numDeriv (optional; for Hessian/Std.Err calculations)

# Minimal example

For example inputs, simulate from a two-population demographic history using `msprime` in python,
writing the simulated tree sequences out to the native binary format of `tskit`:
```{python}
import msprime

TODO
```
Here, we've generated 10 replicate tree sequences (e.g. chromosomes), each 50 Mb and containing genealogies for 20 haploids per population.

In R, first we calculate trio first coalescence rates. We generate 1000 bootstrap replicates by resampling 1000 contiguous blocks of trees, and combine rates across the replicate tree sequences.
```{r}
library(coaldecoder)
TODO
```

Then we need to estimate a precision matrix for the rates from the bootstrap replicates. This is needed to give the model a sense of the natural 'variability' of the rates. In general, if there are a lot of rate statistics (e.g. many populations/epochs) then an enormous number of bootstrap replicates will be needed to get a full-rank precision matrix with low error. So, instead we'll calculate an full-rank approximation. The simplest approximation is a diagonal precision matrix, but this ignores dependence between rate statistics. A better approach is to calculate a shrinkage estimator:
```{r}
corpcor::invcov.shrink(?)
TODO
```

Then we need to set up the model, which requires an array of demographic parameters with dimensions `(# populations, # populations, # epochs)` and the durations of the epochs in generations. The parameterization used is identical to that of `msprime`.
```{r}
TODO
```

Finally, we can fit the model:
```{r}
TODO
```
Here we've specified a smoothing penalty of 1 on the squared differences between log10 parameter values from adjacent epochs. When the temporal resolution is fine, the penalty will have an enormous impact on the fit, and so it's a good idea to use cross-validation to choose the degree of penalization. That could be done using a holdout set of chromosomes, but we skip that here.

`coaldecoder()` returns an object identical in structure to the output of R's `optim()`, so see `?optim` for details. What is relevant here is,
```{r}
str(fitted_model$par) # log10 optimized demographic parameters
str(fitted_model$hessian) # Hessian matrix, if requested
```

To visualise the fitted model, use:
```{r}
TODO
```

Let's that compare to the true model. 
```{r}
true_model <- array(NA, dim=dim(?,?,?))
TODO
```

# Population mergers and identifiability

TODO

# Cross-validation

TODO

# Low-level interface

TODO
