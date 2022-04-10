library(coaldecoder)
library(Matrix)

#---------------no bootstraps, not dividing by time-------------#
P <- 2

y <- array(1, c(12, 3, 1))
n <- matrix(100, 12, 1)
deco <- CoalescentDecoder$new(2, y, n, c(1000.,1000.,1000.), FALSE)

M <- array(0, c(2,2,3))
M[1,1,] <- 1000
M[2,2,] <- 2000
M[1,2,] <- 1e-5
M[2,1,] <- 2e-5

A <- array(0, c(2,2,3))
A[1,1,] <- A[2,2,] <- 0.3
A[1,2,] <- A[2,1,] <- 0.7

X <- deco$initial_state_vectors()
foo <- deco$loglikelihood(X, M, A)

gra_M <- numDeriv::grad(function(M)
               {
               foo <- deco$loglikelihood(X, 10^M, A)
               foo$loglikelihood
}, log10(M))

test_1.1 <- all(dplyr::near(c(foo$gradient$M) * log(10) * c(M), gra_M))

gra_A <- numDeriv::grad(function(A)
               {
               foo <- deco$loglikelihood(X, M, A)
               foo$loglikelihood
}, A)

test_1.2 <- all(dplyr::near(c(foo$gradient$A), gra_A))

gra_X <- numDeriv::grad(function(X)
               {
               foo <- deco$loglikelihood(X, M, A)
               foo$loglikelihood
}, X)

test_1.3 <- all(dplyr::near(c(foo$gradient$X), gra_X))


#------------------____DEPR____-----------------------#

#---------------no bootstraps, not dividing by time-------------#
# testing equivalence with earlier impelementation, so no admixture

y <- array(1, c(12, 3, 1))
n <- matrix(100, 12, 1)
deco <- CoalescentDecoder$new(2, y, n, c(10.,20.,30.), FALSE, TRUE)

A <- array(0, c(2,2,3))
A[1,1,] <- A[2,2,] <- 1
X <- deco$initial_state_vectors()
foo <- deco$loglikelihood(X, M, A)

deco2 <- decoder_gaussian$new(2,3,TRUE)
y2 <- exp(deco2$coalescence_rates(matrix(c(100,200,200,100),4,1),array(1, c(3,4,3))))

test_2.1 <- all(dplyr::near(y2[y2 > 0], deco$observed_rates()[1:6,]))

prec_mat <- diag(12 * 3)
prec_mat[c(y2 == 0),c(y2 == 0)] <- 0
prec_mat <- as(prec_mat, "dgCMatrix")
foo2 <- deco2$deviance(deco2$initial_states(), prec_mat, y2, matrix(0:1, 2, 3), M, c(10,20,30))

test_2.2 <- all(dplyr::near(c(foo$y_hat[1:6,]), foo2$predicted[foo2$predicted > 0]))

res2 <- c(foo2$emission - foo2$predicted)

test_2.3 <- all(dplyr::near(c(foo$residuals[1:6,]), res2[foo2$predicted > 0]))

#likelihoods are not directly comparable cuz there's more states in new model
#as special cast I set it to use pair statistics, then test loglikelihood equivalence via
#test_2.4 <- dplyr::near(foo$loglikelihood*-2, foo2$deviance)

### TODO need more testing.

#test:
# --occupancy function
# --smoothness penalty
# --more populations?
# --divide by time--check gradient

#----------------test with bootstraps, but not dividing rates by time---------------#
boot <- 50
y_boot <- array(1, c(nrow(y),ncol(y),boot+1))
y_boot[,,1] <- y
for(i in 1:boot+1) y_boot[,,i] <- 1 + rnorm(nrow(y), 0, 0.1)
n_boot <- matrix(n, nrow(n), boot+1)
deco_boot <- CoalescentDecoder$new(2, y_boot, n_boot, c(10.,20.,30.), FALSE, FALSE)

M <- array(0, c(2,2,3))
M[1,1,] <- 1000
M[2,2,] <- 2000
M[1,2,] <- 1e-5
M[2,1,] <- 2e-5

A <- array(0, c(2,2,3))
A[1,1,] <- A[2,2,] <- 0.3
A[1,2,] <- A[2,1,] <- 0.7
X <- deco_boot$initial_state_vectors()
foo <- deco_boot$loglikelihood(X, M, A)

#TODO: these are failing but are very close, I think I should just increase tolerance
gra_M <- numDeriv::grad(function(M)
               {
               foo <- deco_boot$loglikelihood(X, 10^M, A)
               foo$loglikelihood
}, log10(M))

test_3.1 <- all(dplyr::near(c(foo$gradient$M) * log(10) * c(M), gra_M))

gra_A <- numDeriv::grad(function(A)
               {
               foo <- deco_boot$loglikelihood(X, M, A)
               foo$loglikelihood
}, A)

test_3.2 <- all(dplyr::near(c(foo$gradient$A), gra_A))

gra_X <- numDeriv::grad(function(X)
               {
               foo <- deco_boot$loglikelihood(X, M, A)
               foo$loglikelihood
}, X)

test_3.3 <- all(dplyr::near(c(foo$gradient$X), gra_X))

#----------------test with bootstraps, dividing rates by time---------------#
#be SURE to check gradient
#TODO:

#----------------test smoothness penalty at different orders---------------#
deco <- CoalescentDecoder$new(2, c(10.,20.,30.,40.,50.), FALSE, FALSE)

M <- array(0, c(2,2,5))
M[1,1,] <- 10^rnorm(dim(M)[3], 3, 0.1)
M[2,2,] <- 10^rnorm(dim(M)[3], 3, 0.1)
M[1,2,] <- 10^rnorm(dim(M)[3], -4, 0.1)
M[2,1,] <- 10^rnorm(dim(M)[3], -4, 0.1)
penalty <- matrix(c(1,0.5,2,0.3),2,2,byrow=TRUE)

foo <- deco$smoothness_penalty(M, penalty, 1)
gra_M <- numDeriv::grad(function(M)
               {
               foo <- deco$smoothness_penalty(10^M, penalty, 1)
               foo$penalty
}, log10(M))

test_5.1 <- all(dplyr::near(c(foo$gradient * log(10) * M), gra_M))

foo <- deco$smoothness_penalty(M, penalty, 2)
gra_M <- numDeriv::grad(function(M)
               {
               foo <- deco$smoothness_penalty(10^M, penalty, 2)
               foo$penalty
}, log10(M))

test_5.2 <- all(dplyr::near(c(foo$gradient * log(10) * M), gra_M))
