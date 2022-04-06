library(coaldecoder)
library(Matrix)

#-------------- test identity admixture matrix --------------#
set.seed(1)
A <- diag(3)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo <- test_TrioAdmixtureProportions(A, G)
sum(foo$matrix) == 63
sum(diag(foo$matrix)) == 63

#-------------- test partial admixture matrix --------------#
set.seed(1)
A <- matrix(c(1,0,0,0.3,0.2,0.5,0,0,1), 3, 3, byrow=TRUE)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo2 <- test_TrioAdmixtureProportions(A, G)

all(dplyr::near(rowSums(foo2$matrix), 1))

#-------------- test full admixture matrix --------------#
set.seed(1)
A <- matrix(c(0.9,0.1,0,0.3,0.2,0.5,0.1,0.1,0.8), 3, 3, byrow=TRUE)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo3 <- test_TrioAdmixtureProportions(A, G)

all(dplyr::near(rowSums(foo2$matrix), 1))

#-------------- test jacobian ---------------------------#
set.seed(1)
A <- matrix(c(0.9,0.05,0.05,0.3,0.2,0.5,0.1,0.1,0.8), 3, 3, byrow=TRUE)
G <- as(matrix(rnorm((4^3 - 1)^2), 4^3 - 1, 4^3 - 1), "dgCMatrix")
foo4 <- test_TrioAdmixtureProportions(A, G)

jac <- numDeriv::jacobian(function(A)
{
  tmp <- test_TrioAdmixtureProportions(A, G)
  c(as.matrix(tmp$matrix))
}, A)

all(dplyr::near(c(foo4$reverse_differentiate), t(jac) %*% c(as.matrix(G))))

