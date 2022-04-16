#fixing matrix epxponential issues

#----------DELETME when resolved

#---more matexp stuff
library(coaldecoder)
library(Matrix)

P <- 3
epoch_duration <- c(1000,1000,1000)
tmp <- CoalescentDecoder$new(3, epoch_duration, TRUE)
S <- length(tmp$emission_states(c("A","B","C")))

y <- matrix(1, S, length(epoch_duration))
B <- lapply(1:ncol(y), function(i) Matrix::Diagonal(nrow(y)))

deco <- CoalescentDecoder$new(3, epoch_duration, TRUE)

M <- array(0, c(P,P,length(epoch_duration)))
M[1,1,] <- 1000
M[2,2,] <- 2000
M[3,3,] <- 3000
M[1,2,] <- 1e-5
M[1,3,] <- 1e-4
M[2,1,] <- 2e-5
M[2,3,] <- 2e-4
M[3,1,] <- 1e-6
M[3,2,] <- 1e-5

A <- array(0, c(P,P,length(epoch_duration)))
set.seed(1); for(i in 1:P) for(j in 1:P) {xx <- rexp(P); A[,i,j] <- xx/sum(xx)}

X <- deco$initial_state_vectors()
X[] <- rexp(length(X[]))
X <- prop.table(X, 2)

foo <- deco$loglikelihood(y, B, X, M, A)

#i would to to figure out how to stabilize this
M2 <- matrix(0, P, P)
M2[upper.tri(M2)] <- 0.1
M2[lower.tri(M2)] <- 1e-8
diag(M2) <- 10
Q <- deco$transition_rates(M2)
yay <- SpMatExp$new(t(Q), foo$X, 9000) 


