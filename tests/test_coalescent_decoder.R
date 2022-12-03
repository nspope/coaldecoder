library(coaldecoder)
library(Matrix)

#---------------likelihood and gradient------------------------#
P <- 3
epoch_duration <- c(1000,1000,1000)
tmp <- CoalescentDecoder$new(3, epoch_duration, TRUE)
S <- length(tmp$emission_states(c("A","B","C")))

y <- matrix(1, S, length(epoch_duration))
B <- matrix(0.5, S, length(epoch_duration))

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

# demographic parameters gradient
gra_M <- numDeriv::grad(function(M)
               {
               foo <- deco$loglikelihood(y, B, X, 10^M, A)
               foo$loglikelihood
}, log10(M))

test_1.1 <- all(dplyr::near(c(foo$gradient$M) * log(10) * c(M), gra_M))

# admixture parameters gradient
arr_softmax <- function(a) {
  for(i in 1:dim(a)[3])
  {
    for(j in 1:dim(a)[2])
    {
      a[,j,i] <- exp(a[,j,i])/sum(exp(a[,j,i]))
    }
  }
  a
}
arr_logit <- function(a) {
  for(i in 1:dim(a)[3])
  {
    for(j in 1:dim(a)[2])
    {
      a[,j,i] <- log(a[,j,i])
    }
  }
  a
}
stopifnot(all(dplyr::near(arr_softmax(arr_logit(A)), A)))

gra_A <- numDeriv::grad(function(A)
               {
               foo <- deco$loglikelihood(y, B, X, M, arr_softmax(A))
               foo$loglikelihood
}, arr_logit(A))

jac_A <- numDeriv::jacobian(arr_softmax, arr_logit(A))
test_1.2 <- all(dplyr::near(t(jac_A) %*% c(foo$gradient$A), gra_A))

# initial states gradient
mat_softmax <- function(a) {
  for(j in 1:dim(a)[2])
  {
    a[,j] <- exp(a[,j])/sum(exp(a[,j]))
  }
  a
}
mat_logit <- function(a) {
  for(j in 1:dim(a)[2])
  {
    a[,j] <- log(a[,j])
  }
  a
}
stopifnot(all(dplyr::near(mat_softmax(mat_logit(X)), X)))

gra_X <- numDeriv::grad(function(X)
               {
               foo <- deco$loglikelihood(y, B, mat_softmax(X), M, A)
               foo$loglikelihood
}, mat_logit(X))

jac_X <- numDeriv::jacobian(mat_softmax, mat_logit(X))
test_1.3 <- all(dplyr::near(t(jac_X) %*% c(foo$gradient$X), gra_X))

stopifnot(all(c(test_1.1, test_1.2, test_1.3)))

#------------------test matrix exponential multiply---------------#
# scaling issues typically happen when: migration rate is very high OR population size is very small

M2 <- matrix(0, P, P)
M2[upper.tri(M2)] <- 0.1
M2[lower.tri(M2)] <- 1e-8
diag(M2) <- 10
tol <- 1e-12

# ensure we have a valid test
# now i've fixed, should work great
X2 <- deco$transition_operator_unsafe(foo$X, M2, 10000)
stopifnot(all(X2 >= 0) & all(X2 <= 1) & all(abs(colSums(X2) - 1) <= tol))

# check against external safe implementation
X3 <- expm::expm(10000 * t(deco$transition_rates(M2))) %*% foo$X
test_2.3 <- all(X3 >= 0) & all(X3 <= 1) & all(abs(colSums(X3) - 1) <= tol)

# compare implementations
test_2.4 <- all(abs(X3 - X2) <= tol)

stopifnot(all(c(test_2.1, test_2.2, test_2.3, test_2.4)))

#---------------------test occupancy function------------------#

#---------------------test penalty function--------------------#
library(coaldecoder)

P <- 3
epoch_duration <- rep(1000, 30)
deco <- CoalescentDecoder$new(P, epoch_duration, TRUE)

set.seed(1)
M <- array(0, c(P,P,length(epoch_duration)))
M[1,1,] <- rnorm(length(epoch_duration), 10000, 100)
M[2,2,] <- rnorm(length(epoch_duration), 10000, 100)
M[3,3,] <- rnorm(length(epoch_duration), 10000, 100)
M[1,2,] <- rnorm(length(epoch_duration), 1e-5, 1e-6)
M[1,3,] <- rnorm(length(epoch_duration), 1e-5, 1e-6)
M[2,1,] <- rnorm(length(epoch_duration), 1e-5, 1e-6)
M[2,3,] <- rnorm(length(epoch_duration), 1e-5, 1e-6)
M[3,1,] <- rnorm(length(epoch_duration), 1e-5, 1e-6)
M[3,2,] <- rnorm(length(epoch_duration), 1e-5, 1e-6)

M[1,1,1:20] <- Inf
M[1,2,1:10] <- 0
M[2,3,10:20] <- 0 #period of isolation
M[3,3,2:30] <- Inf #only present for 1 epoch
M[3,1,1:30] <- 0 #not present

penalty <- matrix(1+rexp(9), 3, 3)

#first order
pen <- deco$smoothness_penalty(M, penalty, 1)
pen_gra <- c(pen$gradient * M * log(10))
pen_gra[is.na(pen_gra)] <- 0
first_ord <- function(x, p){
  z <- x[is.finite(x)]
  sum(-0.5 * diff(z)^2 * p^2)
}
pen_verify <- 0
for(i in 1:3) for(j in 1:3) {
  pen_verify <- pen_verify + first_ord(log10(M[i,j,]), penalty[i,j])
}
gra_M <- numDeriv::grad(function(M) deco$smoothness_penalty(10^M, penalty, 1)$penalty, log10(M))

stopifnot(dplyr::near(pen_verify, pen$penalty))
stopifnot(all(dplyr::near(pen_gra, gra_M)))

#second order
pen <- deco$smoothness_penalty(M, penalty, 2)
pen_gra <- c(pen$gradient * M * log(10))
pen_gra[is.na(pen_gra)] <- 0
second_ord <- function(x, p){
  z <- x[is.finite(x)]
  sum(-0.5 * diff(diff(z))^2 * p^2)
}
pen_verify <- 0
for(i in 1:3) for(j in 1:3) {
  pen_verify <- pen_verify + second_ord(log10(M[i,j,]), penalty[i,j])
}
gra_M <- numDeriv::grad(function(M) deco$smoothness_penalty(10^M, penalty, 2)$penalty, log10(M))

stopifnot(dplyr::near(pen_verify, pen$penalty))
stopifnot(all(dplyr::near(pen_gra, gra_M)))

#---------------------test state labels------------------------#
pop_names <- c("Zebra", "Goat", "Bob")
deco$emission_states(pop_names)
deco$initial_states(pop_names)
deco$transitory_states(pop_names)

#---------------------------------------------------------------#
#TODO: need some more tough edge cases like this one

library(coaldecoder)
library(Matrix)
load("~/Downloads/fing_hell/coaldecoder/problem.RData")

y <- rates
B <- prec
deco <- CoalescentDecoder$new(4, c(t), TRUE)
X <- state
foo <- deco$loglikelihood(y, B, X, M, A)
gra_M <- numDeriv::grad(function(M)
               {
               foo <- deco$loglikelihood(y, B, X, 10^M, A)
               foo$loglikelihood
}, log10(M))

test_1.1 <- all(dplyr::near(c(foo$gradient$M) * log(10) * c(M), gra_M)) #fails, this is really close though

agh <- deco$transition_operator(X, M[,,1], t)
erg <- expm::expm(t * t(deco$transition_rates(M[,,1]))) %*% X
test_1.2 <- all(dplyr::near(agh, erg))

#----------------keep this as a test case---------------#
#TODO need a suite of tough test cases
library(coaldecoder)
library(Matrix)
load("~/Downloads/fing_hell/problem_dem.RData")
deco <- CoalescentDecoder$new(4, c(t), TRUE)
X <- state
Q <- t(deco$transition_rates(M[,,1]))
Spoo <- SpMatExp$new(Q, X, t)
stopifnot(all(dplyr::near( Spoo$result, (expm(t * Q) %*% X) )))

#---------------test one-norm estimation-----------------#
#TODO have saved test cases
#also test exact norm on small example
library(coaldecoder)
library(Matrix)

set.seed(1)
p <- 4
eek <- rsparsematrix(10000, 10000, 0.01)
Ax = function(x) {for(i in 1:p) x <- eek %*% x; x}
Atx = function(x) {for(i in 1:p) x <- t(eek) %*% x; x}
ah <- onenormest(A.x=Ax, At.x=Atx, n=nrow(eek))
ah$est

writeMM(eek, "~/Dropbox/Projects/expm/rsparsematrix.mm")

set.seed(1)
eh <- OneNormEst$new(eek, 4, TRUE)
eh$est
