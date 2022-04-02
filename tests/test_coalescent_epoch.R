library(coaldecoder)
library(Matrix)

# --------- test equivalence with pretested mockup ---------#
set.seed(1)
A <- diag(3)
M <- matrix(c(1e4,1e-4,1e-5,1e-6,5e4,1e-3,2e-4,3e-4,5e4),3,3,byrow=TRUE)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo <- test_TrioTransitionRates(M, G)
states <- matrix(0, nrow(G), ncol(foo$i_map_s))
for(i in 1:ncol(foo$i_map_s)) states[foo$i_map_s[2,i]+1, foo$i_map_s[1,i]+1] <- 1

gr <- matrix(rnorm(prod(dim(states))), nrow(states), ncol(states))
foo <- test_CoalescentEpoch(states, M, A, gr)

deco <- decoder_gaussian$new(3,3,TRUE)
states2 <- deco$initial_states()
gr2 <- matrix(rnorm(prod(dim(states2))), nrow(states2), ncol(states2))
bar <- test_coalescent_epoch_gaussian(states2, M, gr2)
y_hat2 <- t(bar$y_hat)
y_hat2 <- c(y_hat2[y_hat2 > 0])
resid2 <- t(bar$residual)
resid2 <- resid2[resid2 > 0]

test_1.1 <- all(dplyr::near(foo$y_hat[1:18], y_hat2))
test_1.2 <- all(dplyr::near(foo$residual[1:18], resid2))

# --------- test gradient w/o admixture ---------#
set.seed(1)
A <- diag(3)
M <- matrix(c(1e4,1e-4,1e-5,1e-6,5e4,1e-3,2e-4,3e-4,5e4),3,3,byrow=TRUE)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo <- test_TrioTransitionRates(M, G)
states <- matrix(0, nrow(G), ncol(foo$i_map_s))
for(i in 1:ncol(foo$i_map_s)) states[foo$i_map_s[2,i]+1, foo$i_map_s[1,i]+1] <- 1

gr <- matrix(rnorm(prod(dim(states))), nrow(states), ncol(states))
foo2 <- test_CoalescentEpoch(states, M, A, gr)

#gradient wrt input states
gra_states <- numDeriv::grad(function(states){
  foo2 <- test_CoalescentEpoch(states, M, A, gr)
 foo2$loglikelihood
}, states)
jac_states <- numDeriv::jacobian(function(states){
  foo2 <- test_CoalescentEpoch(states, M, A, gr)
  foo2$states_end
}, states)

test_2.1 <- all(dplyr::near(c(foo2$gradient_states), t(jac_states) %*% c(gr) + gra_states))

#gradient wrt migration matrix
gra_M <- numDeriv::grad(function(M){
  foo2 <- test_CoalescentEpoch(states, 10^M, A, gr)
  foo2$loglikelihood
}, log10(M))
jac_M <- numDeriv::jacobian(function(M){
  foo2 <- test_CoalescentEpoch(states, 10^M, A, gr)
  foo2$states_end
}, log10(M))

test_2.2 <- all(dplyr::near(c(foo2$gradient_M) * c(M) * log(10), t(jac_M) %*% c(gr) + gra_M))

# --------- test w/ admixture (w/ gradient) -----#
set.seed(1)
A <- matrix(c(0.1,0.6,0.3,0.2,0.4,0.4,0.8,0.1,0.1),3,3,byrow=3)
M <- matrix(c(1e4,1e-4,1e-5,1e-6,5e4,1e-3,2e-4,3e-4,5e4),3,3,byrow=TRUE)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo <- test_TrioTransitionRates(M, G)
states <- matrix(0, nrow(G), ncol(foo$i_map_s))
for(i in 1:ncol(foo$i_map_s)) states[foo$i_map_s[2,i]+1, foo$i_map_s[1,i]+1] <- 1

gr <- matrix(rnorm(prod(dim(states))), nrow(states), ncol(states))
foo3 <- test_CoalescentEpoch(states, M, A, gr)

#gradient wrt input states
gra_states <- numDeriv::grad(function(states){
  foo3 <- test_CoalescentEpoch(states, M, A, gr)
 foo3$loglikelihood
}, states)
jac_states <- numDeriv::jacobian(function(states){
  foo3 <- test_CoalescentEpoch(states, M, A, gr)
  foo3$states_end
}, states)

test_3.1 <- all(dplyr::near(c(foo3$gradient_states), t(jac_states) %*% c(gr) + gra_states))

#gradient wrt migration matrix
gra_M <- numDeriv::grad(function(M){
  foo3 <- test_CoalescentEpoch(states, 10^M, A, gr)
  foo3$loglikelihood
}, log10(M))
jac_M <- numDeriv::jacobian(function(M){
  foo3 <- test_CoalescentEpoch(states, 10^M, A, gr)
  foo3$states_end
}, log10(M))

test_3.2 <- all(dplyr::near(c(foo3$gradient_M) * c(M) * log(10), t(jac_M) %*% c(gr) + gra_M))

#gradient wrt admixture matrix
gra_A <- numDeriv::grad(function(A){
  foo3 <- test_CoalescentEpoch(states, M, 10^A, gr)
  foo3$loglikelihood
}, log10(A))
jac_A <- numDeriv::jacobian(function(A){
  foo3 <- test_CoalescentEpoch(states, M, 10^A, gr)
  foo3$states_end
}, log10(A))
test_3.3 <- all(dplyr::near(c(foo3$gradient_A) * c(A) * log(10), t(jac_A) %*% c(gr) + gra_A))

#-----------------
stopifnot(all(c(test_1.1, test_1.2, test_1.3, test_2.1, test_2.2, test_3.1, test_3.2, test_3.3)))
#-----------------
