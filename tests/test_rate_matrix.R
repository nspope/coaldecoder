library(coaldecoder)
library(Matrix)

# --------- check valid rate matrix ---------
M <- matrix(1, 3, 3)
set.seed(1); G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo <- test_TrioTransitionRates(M, G)

all(dplyr::near(rowSums(foo$matrix), 0))

# --------- check that gives the same values as old implementation ----------
M <- matrix(c(1e4,1e-4,1e-5,1e-6,5e4,1e-3,2e-4,3e-4,5e4),3,3,byrow=TRUE)
set.seed(1); G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
new_i_state <- '{0,0,1}'; old_i_state <- c(2,1,0)
new_e_state <- c('t1::((0,1),0)','t2::((0,1),0)'); old_e_state <- c(0,1)

foo <- test_TrioTransitionRates(M, G)
mapper <- foo$s_map_e[,foo$s_map_e[1,]==(which(foo$initial_states==new_i_state)-1)]
mapper <- mapper[,mapper[3,] %in% (which(foo$emission_states %in% new_e_state)-1)]
eye <- diag(nrow(foo$matrix))[,which(foo$initial_states==new_i_state)]
a <- expm::expm(t(foo$matrix)) %*% eye
new_prob <- sum(a[mapper[2,]+1])

deco <- decoder_gaussian$new(3,3,TRUE)
bar <- deco$transition_probabilities(c(0,1,2), M, 1)
old_i_state_index <- which(apply(deco$emission_classes(), 1, paste, collapse=",") == paste(old_i_state, collapse=","))
eye2 <- diag(nrow(bar))[,old_i_state_index]
b <- bar %*% eye2
last3 <- (length(eye2):(length(eye2)-2))[3:1]
c_map <- deco$map_coalescent_states_to_emission_states()[,old_i_state_index] + 1
old_e_state_index <- which(apply(deco$emission_states()[c_map,], 1, paste, collapse=",") == paste(old_e_state, collapse=","))
old_prob <- sum(b[last3[old_e_state_index]])

dplyr::near(new_prob, old_prob)

#------------ check reverse differentiation ----------
M <- matrix(c(1e4,1e-4,1e-5,1e-6,5e4,1e-3,2e-4,3e-4,5e4),3,3,byrow=TRUE)
set.seed(1); G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
jac <- numDeriv::jacobian(function(M){
  foo <- test_TrioTransitionRates(M, G)
  c(as.matrix(foo$matrix))
}, M)
foo <- test_TrioTransitionRates(M, G)

all(dplyr::near(c(foo$reverse_differentiate), t(jac) %*% c(as.matrix(G))))

#------------- check rate matrix for "ancient" populations ------------
#for this, we set effective population size to Inf and migration to 0 this
#should keep state probabilities for all states involving the "ancient"
#population constant when multiplied by the transition matrix

# two pops
M <- matrix(c(1e4,0,0,Inf),2,2,byrow=TRUE)
deco <- CoalescentDecoder$new(nrow(M), c(1000), TRUE, TRUE)
start <- deco$initial_state_vectors()
end <- expm::expm(1000*t(deco$transition_rates(M))) %*% start
const <- sapply(1:ncol(start), function(i) all(start[,i] == end[,i]))

#because lineages in first pop can't migrate, state vectors where pop1 can't coalesce are constant
should_be_const <- c("{0,1,1}", "{1,1,1}")
test_4.1 <-  all(const[deco$initial_states() %in% should_be_const]) & sum(const)==2

# three pops
M <- matrix(c(1e4,1e-4,0,1e-6,5e4,0,0,0,Inf),3,3,byrow=TRUE)
deco <- CoalescentDecoder$new(nrow(M), c(1000), TRUE, TRUE)
start <- deco$initial_state_vectors()
end <- expm::expm(1000*t(deco$transition_rates(M))) %*% start
const <- sapply(1:ncol(start), function(i) all(start[,i] == end[,i]))

#because lineages in first pop CAN migrate, only state vector with all samples in pop 3 is constant
should_be_const <- c('{2,2,2}')
test_4.2 <- all(const[deco$initial_states() %in% should_be_const]) & sum(const)==1

