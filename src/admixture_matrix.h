//
// MIT License
//
// Copyright (c) 2021-2022 Nathaniel S. Pope
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#ifndef _ADMIXTURE_MATRIX
#define _ADMIXTURE_MATRIX

#include <RcppArmadillo.h> 
#include <vector>
#include <string>
#include <algorithm>

// 1. analyze structure of admixture matrix to get "contributors" to each population
//    for example we if we have 3,4 -> 0
//    then we need a data structure that gives x[0] = {3, 4}
// 2. loop over all states
// 3. loop over lineages in states
// 4. visit populations in the data structure based on current location of lineage.
//    so for example, lets say we are at state u = {0, 0, 4}
//    and we have x[0] = {3,4} and x[4] = {}
//    then, we visit states v = {a, b, 4} where a,b %in% {3,4}
// 5. the "contribution" weight for v is product of admixture proportions, that is
//    A[0,a] * A[0,b]

struct TrioAdmixtureProportions
{
  /*
   *  Transforms trio state probabilities according to population admixture
   *  matrix A, where A[i,j] is the proportion of population i that comes from
   *  population j.
   *
   *  This remains relatively computationally cheap even when A as large, as
   *  long as most entries are 0 (e.g. A is sparse).
   */

  const std::string prefix = "[TrioAdmixtureProportions] ";
  const bool check_valid = true;

  const unsigned P; //number of populations
  const unsigned S; //number of states
  const arma::mat A; //admixture parameter matrix
  const arma::sp_mat X; //transition rate matrix

  TrioAdmixtureProportions (const arma::mat& _A, const bool _check_valid = true)
    : check_valid (_check_valid)
    , P (_A.n_rows)
    , S (power(P+1, 3) - 1)
    , A (_A)
    , X (build_matrix())
  {}

  unsigned power(unsigned x, unsigned p) const
  {
    if (p == 0) return 1;
    if (p == 1) return x;
          
    unsigned tmp = power(x, p/2);
    if (p%2 == 0) return tmp * tmp;
    else return x * tmp * tmp;
  }

  arma::uword linear_index (const arma::uvec::fixed<3>& populations) const
  {
    if (arma::any(populations > P))
    {
      Rcpp::stop(prefix + " Linear index out-of-range");
    }

    arma::uword index = 0;
    for(int i=2; i>=0; --i)
    {
      arma::uword partial_index = populations[i];
      for(int j=2; j>i; --j)
      {
        partial_index *= (P + 1);
      }
      index += partial_index;
    }
    return index;
  }

  void push_back (std::vector<arma::uword>& R, std::vector<arma::uword>& C, std::vector<double>& V, const arma::uword& r, const arma::uword& c, const double& v) const
  {
    R.push_back(r);
    C.push_back(c);
    V.push_back(v);
  }

  std::vector<arma::uvec> sparsity_pattern (void) const
  {
    /*
     *  Analyze sparsity of admixture matrix
     */

    if (A.n_rows != P || A.n_cols != P)
    {
      Rcpp::stop(prefix + "Admixture matrix dimensions must equal populations");
    }

    if (check_valid)
    {
      if (arma::any(arma::vectorise(A) < 0.0) || arma::any(arma::vectorise(A) > 1.0))
      {
        Rcpp::stop(prefix + "Admixture proportions invalid");
      }
    }

    std::vector<arma::uvec> contributors (P + 1);
    for (unsigned i=0; i<P; ++i)
    {
      if (check_valid)
      {
        if (std::fabs(arma::accu(A.col(i)) - 1.0) > A.n_rows*arma::datum::eps)
        {
          Rcpp::stop(prefix + "Admixture proportions do not sum to one");
        }
      }
      contributors[i] = arma::find(A.row(i) > 0.0);
    }
    contributors[P] = arma::uvec({P});

    return contributors;
  }

  arma::sp_mat build_matrix (void) const
  {
    std::vector<arma::uvec> contributors = sparsity_pattern();

    arma::uvec::fixed<3> u;

    // number of non-zeros
    arma::uword storage_bound = 0;
    for (u[0]=0; u[0]<contributors.size(); ++u[0])
    {
      for (u[1]=0; u[1]<contributors.size(); ++u[1])
      {
        for (u[2]=0; u[2]<contributors.size(); ++u[2])
        {
          for (auto i : contributors[u[0]])
          {
            for (auto j : contributors[u[1]])
            {
              for (auto k : contributors[u[2]])
              {
                storage_bound += 1;
              }
            }
          }
        }
      }
    }
    storage_bound -= 1;

    std::vector<arma::uword> cols;
    std::vector<arma::uword> rows;
    std::vector<double> values;
    cols.reserve(storage_bound);
    rows.reserve(storage_bound);
    values.reserve(storage_bound);

    for (u[0]=0; u[0]<contributors.size(); ++u[0])
    {
      for (u[1]=0; u[1]<contributors.size(); ++u[1])
      {
        for (u[2]=0; u[2]<contributors.size(); ++u[2])
        {
          if (!arma::all(u == P))
          {
            arma::uvec::fixed<3> _u = u;
            arma::uword dest_index = linear_index(u);

            for (auto i : contributors[u[0]])
            {
              for (auto j : contributors[u[1]])
              {
                for (auto k : contributors[u[2]])
                {
                  _u[0] = i; _u[1] = j; _u[2] = k;
                  arma::uword source_index = linear_index(_u);

                  double rate = 1.0;
                  for (unsigned l=0; l<3; ++l)
                  {
                    if (u[l] < P && _u[l] < P)
                    {
                      rate *= A(u[l],_u[l]);
                    }
                  }

                  push_back(rows, cols, values, dest_index, source_index, rate);
                }
              }
            }
          }
        }
      }
    }

    if (values.size() != storage_bound)
    {
      std::cout << storage_bound << " " << values.size() << std::endl;
      Rcpp::stop(prefix + " Storage bound violated");
    }

    arma::umat locations = arma::join_vert(
        arma::urowvec(rows),
        arma::urowvec(cols)
    );
    arma::sp_mat out (true, locations, arma::vec(values), S, S);
    return out;
  }  
  
  arma::mat reverse_differentiate (const arma::sp_mat& gradient) const
  {
    std::vector<arma::uvec> contributors = sparsity_pattern();

    arma::mat out = arma::zeros(P, P);

    arma::uvec::fixed<3> u;
    for (u[0]=0; u[0]<contributors.size(); ++u[0])
    {
      for (u[1]=0; u[1]<contributors.size(); ++u[1])
      {
        for (u[2]=0; u[2]<contributors.size(); ++u[2])
        {
          if (!arma::all(u == P))
          {
            arma::uvec::fixed<3> _u = u;
            arma::uword dest_index = linear_index(u);

            for (auto i : contributors[u[0]])
            {
              for (auto j : contributors[u[1]])
              {
                for (auto k : contributors[u[2]])
                {
                  _u[0] = i; _u[1] = j; _u[2] = k;
                  arma::uword source_index = linear_index(_u);

                  for (unsigned l=0; l<3; ++l)
                  {
                    if (u[l] < P && _u[l] < P)
                    {
                      out.at(u[l],_u[l]) += 
                        gradient.at(dest_index, source_index) *
                        X.at(dest_index, source_index) / A.at(u[l],_u[l]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    return out;
  }
};

#endif
