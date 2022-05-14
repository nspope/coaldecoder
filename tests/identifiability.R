##---------- identifiability --------------#

##a "data agnostic" measure of identifiability,
## \sum_{topology} (\partial rates_{topology} / \partial parameter)^2
##e.g. norm of jacobian with respect to a given parameter
#rate_jacobian <- numDeriv::jacobian(function(M) 
#{
#  population_tree$set_demographic_parameters(10^M)
#  log10(population_tree$expected_coalescence_rates())
#}, log10(true_parameters_coarse))
#
#identifiability <- array(apply(rate_jacobian, 2, function(x) sum(x^2)), dim(true_parameters_coarse))
