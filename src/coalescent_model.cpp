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
#include <RcppArmadillo.h> 
#include <vector>
#include <string>
#include "expm_revdiff_sparse.h"
#include "transition_rate_matrix.h"
#include "admixture_matrix.h"

struct CoalescentEpoch
{
  const std::string prefix = "[CoalescentEpoch] ";
  const bool check_valid = true;
  const bool use_marginal_statistics = true; //see note in cnstr

  //TODO this should depend on dimension of state vectors
  const double left_stochastic_tol = 1e-12; 

  const TrioAdmixtureProportions admix_matrix;
  const TrioTransitionRates rate_matrix;

  const arma::umat emission_mapping;
  const arma::umat state_mapping;
  const arma::umat initial_mapping;

  const arma::vec y; //observed statistics
  const arma::vec B; //bootstrap sqrt(precision) of y
  const arma::mat A; //admixture proportions
  const arma::mat M; //demographic parameters
  const double t; //duration of epoch

  arma::vec y_hat;
  arma::vec residual;
  arma::vec gradient;
  arma::vec uncoalesced; //uncoalesced lineages at start
  arma::mat states; //states immediately after admixture
  double loglikelihood;

  CoalescentEpoch (
      arma::mat& _states,
      const arma::vec& _y, 
      const arma::vec& _B,
      const arma::mat& _M,
      const arma::mat& _A,
      const double& _t,
      // mappings
      const arma::sp_mat& _rate_matrix_template,
      const arma::umat& _emission_mapping,
      const arma::umat& _state_mapping,
      const arma::umat& _initial_mapping,
      const bool _check_valid = true
  ) : check_valid (_check_valid)
    , admix_matrix (_A, _check_valid)
    , rate_matrix (_M, _rate_matrix_template, _check_valid)
    , emission_mapping (_emission_mapping)
    , state_mapping (_state_mapping)
    , initial_mapping (_initial_mapping)
    , y (_y)
    , B (_B)
    , A (_A)
    , M (_M)
    , t (_t)
  {
    /* 
     *  Update states and calculate likelihood, gradient of emissions within epoch.
     */

    if (_y.n_elem != emission_mapping.n_cols)
    {
      Rcpp::stop(prefix + "Vector of rates has the wrong dimension");
    }

    if (_B.n_elem != emission_mapping.n_cols )
    {
      Rcpp::stop(prefix + "Bootstrap precision has the wrong dimension");
    }

    if (_states.n_rows != arma::accu(rate_matrix.S) || 
        _states.n_cols != initial_mapping.n_cols)
    {
      Rcpp::stop(prefix + "State probability vectors have the wrong dimension");
    }

    if (check_valid)
    {
      if (_t <= 0.0)
      {
        Rcpp::stop(prefix + "Epoch duration must be positive");
      }

      if (arma::any(arma::abs(arma::sum(_states, 0) - 1.0) > left_stochastic_tol))
      {
        Rcpp::stop(prefix + "State probability vector does not sum to one");
      }
      if (arma::any(arma::vectorise(_states) < 0.0))
      {
        Rcpp::stop(prefix + "Negative state probabilities");
      }
    }

    unsigned num_emission = emission_mapping.n_cols;
    unsigned num_initial = initial_mapping.n_cols;

    // Admixture
    states = _states;
    _states = admix_matrix.X * _states;

    // Transition probabilities
    // TODO: this is *still* failing without safeguards
    //SparseMatrixExponentialMultiply transition (arma::trans(rate_matrix.X), _states, t);
    SparseMatrixExponentialMultiplySafe transition (arma::trans(rate_matrix.X), _states, t);
    _states = transition.result;

    // Fitted rates
    arma::vec transitory_start = arma::ones(num_initial);
    arma::vec coalesced_start = arma::zeros(num_emission);
    arma::vec coalesced_end = arma::zeros(num_emission);
    for (unsigned i=0; i<state_mapping.n_cols; ++i)
    {
      arma::uword col = state_mapping.at(0,i);
      arma::uword row = state_mapping.at(1,i);
      arma::uword j = state_mapping.at(2,i);
      transitory_start.at(col) += -states.at(row, col);
      coalesced_start.at(j) += states.at(row, col);
      coalesced_end.at(j) += _states.at(row, col);
    }

    arma::vec subtransitory_start = arma::zeros(num_initial);
    uncoalesced = arma::zeros(num_emission);
    int offset = num_emission / 2;
    for (int j=0; j<num_emission; ++j)
    {
      arma::uword col = emission_mapping.at(0,j);
      arma::uword lin = emission_mapping.at(1,j);

      // sum mass over transitory 2-lineage states in a column
      subtransitory_start[col] += lin == 2 ? coalesced_start[j] : 0.0;

      // 2-lineage states can transition to 1-lineage states
      uncoalesced[j] = lin == 2 ? transitory_start[col] :
        transitory_start[col] + subtransitory_start[col]; 

      // if the input probabilities are marginal, add 1-lineage to 2-lineage.
      // this is because the state vector contains:
      //    p(t1 < x & t2 > x) ===> 2-lineage statistics
      //    p(t1 < x & t2 < x) ===> 1-lineage statistics
      // and we want:
      //    p(t1 < x) = p(t1 < x & t2 < x) + p(t1 < x & t2 > x)
      // which will ensure that the statistics are strictly positive.
      if (lin == 1 && use_marginal_statistics)
      {
        coalesced_start[j - offset] += coalesced_start[j];
        coalesced_end[j - offset] += coalesced_end[j];
      }
    }

    y_hat = (coalesced_end - coalesced_start) / uncoalesced;
    y_hat /= t;

    // Calculate loglikelihood
    residual = B % (y - y_hat);
    gradient = B % residual;
    loglikelihood = -0.5 * arma::dot(residual, residual);
  }

  arma::cube reverse_differentiate (arma::mat& _states)
  {
    if (_states.n_rows != arma::accu(rate_matrix.S) || 
        _states.n_cols != initial_mapping.n_cols)
    {
      Rcpp::stop(prefix + "State gradient vectors have the wrong dimension");
    }

    unsigned num_emission = emission_mapping.n_cols;
    unsigned num_initial = initial_mapping.n_cols;

    // Gradient wrt fitted rates
    arma::vec d_y_hat = gradient * 1.0/t;

    // Gradient wrt state vectors
    arma::mat d_states = arma::zeros(arma::size(states));
    arma::vec d_uncoalesced = -t * y_hat % d_y_hat/uncoalesced;
    arma::vec d_coalesced_start = -d_y_hat/uncoalesced;
    arma::vec d_coalesced_end = d_y_hat/uncoalesced;
    arma::vec d_transitory_start = arma::zeros(num_initial);
    arma::vec d_subtransitory_start = arma::zeros(num_initial);
    int offset = num_emission / 2;
    for (int j=num_emission-1; j>=0; --j)
    {
      arma::uword col = emission_mapping.at(0,j);
      arma::uword lin = emission_mapping.at(1,j);
      if (lin == 1 && use_marginal_statistics)
      {
        d_coalesced_start[j] += d_coalesced_start[j - offset];
        d_coalesced_end[j] += d_coalesced_end[j - offset];
      }
      d_transitory_start[col] += d_uncoalesced[j];
      if (lin == 1) 
      {
        //d_coalesced_start[j - offset] += d_uncoalesced[j];
        d_subtransitory_start[col] += d_uncoalesced[j];
      } else if (lin == 2) {
        d_coalesced_start[j] += d_subtransitory_start[col];
      }
    }
    for (unsigned i=0; i<state_mapping.n_cols; ++i)
    {
      arma::uword col = state_mapping.at(0,i);
      arma::uword row = state_mapping.at(1,i);
      arma::uword j = state_mapping.at(2,i);
      d_states.at(row, col) += -d_transitory_start.at(col);
      d_states.at(row, col) += d_coalesced_start.at(j);
      _states.at(row, col) += d_coalesced_end.at(j);
    }

    // Gradient wrt starting state vectors, trio rate matrix
    // (this adds gradient contribution to existing d_states)
    SparseMatrixExponentialMultiplySafe transition (arma::trans(rate_matrix.X), admix_matrix.X * states, t);
    //SparseMatrixExponentialMultiply transition (arma::trans(rate_matrix.X), admix_matrix.X * states, t);
    arma::mat tmp;
    arma::sp_mat d_rate_matrix = 
      transition.reverse_differentiate(tmp, _states); 
    _states = tmp;

    // Gradient wrt pre-admixture state vectors, trio admixture matrix
    arma::sp_mat d_admix_matrix (arma::size(admix_matrix.X));
    for (arma::sp_mat::const_iterator it = admix_matrix.X.begin(); 
         it != admix_matrix.X.end(); ++it)
    {
      double val = arma::dot(_states.row(it.row()), states.row(it.col()));
      d_admix_matrix.at(it.row(), it.col()) = val;
    }
    _states = arma::trans(admix_matrix.X) * _states;
    _states = _states + d_states;

    // Gradient wrt demographic parameter, admixture matrices
    arma::cube d_parameters (rate_matrix.P, rate_matrix.P, 2);
    d_parameters.slice(0) = rate_matrix.reverse_differentiate(arma::trans(d_rate_matrix));
    d_parameters.slice(1) = admix_matrix.reverse_differentiate(d_admix_matrix);

    return d_parameters;
  }
};

struct CoalescentDecoder
{
  const std::string prefix = "[CoalescentDecoder] ";
  const bool check_valid = true;
  const bool use_marginal_statistics = true; //see notes in method coalescence_rate

  const unsigned P; //number of populations
  const unsigned T; //number of epochs
  const arma::vec t; // epoch duration
  const TrioTransitionRates rate_matrix_template;
  const arma::umat emission_mapping;
  const arma::umat state_mapping;
  const arma::umat initial_mapping;

  CoalescentDecoder (
      const unsigned& _P, // number of populations
      const arma::vec& _t, // epoch duration
      const bool _check_valid = true
  ) : check_valid (_check_valid)
    , P (_P)
    , T (_t.n_elem)
    , t (_t)
    , rate_matrix_template (arma::ones(_P, _P))
    , emission_mapping (rate_matrix_template.emission_to_initial())
    , state_mapping (rate_matrix_template.states_to_emission())
    , initial_mapping (rate_matrix_template.initial_to_states())
  {}

  arma::mat coalescence_rates (const arma::mat& _y, const arma::vec& _n, const arma::vec& _t, const bool _use_rates = true, const bool _use_marginal_statistics = true)
  {
    if (_y.n_rows != emission_mapping.n_cols)
    {
      Rcpp::stop(prefix + "Need a coalescence count for each emission");
    }
    if (_n.n_elem != _y.n_rows || _y.n_cols != _t.n_elem)
    {
      Rcpp::stop(prefix + "Coalescence count dimensions don't match");
    }

    // get denominator; e.g. total mass associated with starting
    // configuration
    arma::mat n = arma::zeros (initial_mapping.n_cols, 2);
    for (unsigned i=0; i<emission_mapping.n_cols; ++i)
    {
      arma::uword col = emission_mapping.at(0,i);
      arma::uword lin = emission_mapping.at(1,i);
      n.at(col, lin-1) += _n.at(i);
    }

    // calculate rates
    arma::mat y (arma::size(_y));
    for (unsigned i=0; i<T; ++i)
    {
      // if not using marginal rates, things get more complicated
      // so let's hold off on that for the moment
      for (unsigned j=0; j<emission_mapping.n_cols; ++j)
      {
        arma::uword col = emission_mapping.at(0,j);
        arma::uword lin = emission_mapping.at(1,j);
        y.at(j, i) = _y.at(j, i) / n.at(col, lin-1);
      }
      for (unsigned j=0; j<emission_mapping.n_cols; ++j)
      {
        arma::uword col = emission_mapping.at(0,j);
        arma::uword lin = emission_mapping.at(1,j);
        n.at(col, lin-1) -= _y.at(j, i);
      }
      if (_use_rates) y.col(i) /= _t.at(i);
    }

    return y;
  }

  std::vector<std::string> initial_states (const std::vector<std::string>& _names) const
  {
    return rate_matrix_template.initial_states(_names);
  }

  std::vector<std::string> emission_states (const std::vector<std::string>& _names) const
  {
    return rate_matrix_template.emission_states(_names);
  }

  std::vector<std::string> transitory_states (const std::vector<std::string>& _names) const
  {
    return rate_matrix_template.transitory_states(_names);
  }

  arma::sp_mat transition_rates (const arma::mat& _M)
  {
    /*
     *  Return trio transition rate matrix given demographic parameters
     */

    TrioTransitionRates rates (_M, rate_matrix_template.X, check_valid);
    return rates.X;
  }

  arma::sp_mat admixture_proportions (const arma::mat& _A)
  {
    /*
     *  Return trio admixture proportions given population admixture
     *  proportions
     */

    TrioAdmixtureProportions admix (_A, check_valid);
    return admix.X;
  }

  arma::mat transition_operator_unsafe (const mat& _X, const arma::mat& _M, const double& _t)
  {
    /*
     *  Map state probability vectors across duration, given demographic
     *  parameters
     */
    TrioTransitionRates rates (_M, check_valid);
    SparseMatrixExponentialMultiply expm (arma::trans(rates.X), _X, _t);
    return expm.result;
  }

  arma::mat transition_operator (const mat& _X, const arma::mat& _M, const double& _t)
  {
    /*
     *  Map state probability vectors across duration, given demographic
     *  parameters
     */
    TrioTransitionRates rates (_M, check_valid);
    SparseMatrixExponentialMultiplySafe expm (arma::trans(rates.X), _X, _t);
    return expm.result;
  }

  arma::mat admixture_operator (const mat& _X, const arma::mat& _A)
  {
    /*
     *  Map state probability vectors given admixture proportions
     */

    TrioAdmixtureProportions admix (_A, check_valid);
    return admix.X * _X;
  }

  arma::cube occupancy_probabilities (const arma::cube& _M, const arma::cube& _A)
  {
    /*
     *  Computes the probability that a lineage starting in a given population
     *  is in other populations at the beginning of epochs
     */

    if (_M.n_slices != T || _M.n_rows != P || _M.n_cols != P)
    {
      Rcpp::stop(prefix + "Demographic parameter array has wrong dimensions");
    }
    if (_A.n_slices != T || _A.n_rows != P || _A.n_cols != P)
    {
      Rcpp::stop(prefix + "Demographic parameter array has wrong dimensions");
    }
    if (check_valid)
    {
      if (arma::any(arma::vectorise(_M) < 0.0))
      {
        Rcpp::stop(prefix + "Demographic parameter array has negative values");
      }
      if (arma::any(arma::vectorise(_A) < 0.0))
      {
        Rcpp::stop(prefix + "Admixture proportions have negative values");
      }
      for (unsigned i=0; i<P; ++i)
      {
        for (unsigned j=0; j<T; j++)
        {
          if (std::fabs(arma::accu(_A.slice(j).col(i)) - 1.0) > _A.n_rows*arma::datum::eps)
          {
            Rcpp::stop(prefix + "Admixture proportions do not sum to one");
          }
        }
      }
    }

    arma::cube out = arma::zeros<arma::cube>(P, P, T + 1);
    out.slice(0).eye();
    for (unsigned i=0; i<T; ++i)
    {
      arma::mat rate_matrix = _M.slice(i);
      for (unsigned j=0; j<P; ++j)
      {
        rate_matrix.at(j,j) = 0.;
        rate_matrix.at(j,j) = -arma::accu(rate_matrix.row(j));
      }
      out.slice(i+1) = 
        arma::expmat(t.at(i) * arma::trans(rate_matrix)) * _A.slice(i) * out.slice(i);
    }

    return out;
  }

  Rcpp::List smoothness_penalty (const arma::cube& _M, const arma::mat& _penalty, const unsigned _order = 1)
  {
    /*
     *  Calculates a penalty on squared differences of a given order, for each
     *  demographic parameter. The differencing is done with respect to time
     *  and on a log10 scale.  Also calculates the gradient wrt the penalty.
     *  Parameters that are are set to Inf or 0 are skipped. 
     *
     *  TODO: if there are splits and mergers involving the same populations,
     *  then there will be "gaps" in the parameter trajectories. These gaps are
     *  ignored; difference between parameter values at the beginning/end of gap
     *  will be penalized!
     */

    if (_order == 0)
    {
      Rcpp::stop(prefix + "Order of differencing operator must be a positive integer");
    }
    if (_M.n_slices != T || _M.n_rows != P || _M.n_cols != P)
    {
      Rcpp::stop(prefix + "Demographic parameter array has wrong dimensions");
    }
    if (arma::any(arma::vectorise(_M) < 0.0))
    {
      Rcpp::stop(prefix + "Demographic parameter array has negative values");
    }
    if (_penalty.n_rows != P || _penalty.n_cols != P)
    {
      Rcpp::stop(prefix + "Penalty matrix has wrong dimensions");
    }
    if (arma::any(arma::vectorise(_penalty) < 0.0))
    {
      Rcpp::stop(prefix + "Penalty matrix has negative values");
    }

    arma::cube gradient_M (arma::size(_M), arma::fill::zeros);

    // OLD using all values
    //// construct difference operator
    //arma::sp_mat diff = arma::speye(T, T);
    //for (unsigned i=0; i<_order; ++i)
    //{
    //  arma::sp_mat d = -1.0 * arma::speye(diff.n_rows - 1, diff.n_rows);
    //  for(unsigned j=0; j<d.n_rows; ++j)
    //  {
    //    d.at(j,j+1) = 1.0;
    //  }
    //  diff = d * diff;
    //}

    // prior, gradient of prior wrt parameters
    double penalty = 0;
    for (unsigned i=0; i<P; ++i)
    {
      for (unsigned j=0; j<P; ++j)
      {
        arma::vec par_ij = arma::vectorise(_M.tube(i,j));
        arma::vec gra_ij = arma::zeros(T);

        // find valid parameter values
        arma::uvec valid;
        if (i == j) 
        { 
          valid = arma::find_finite(par_ij);
        } else {
          valid = arma::find(par_ij > 0.0);
        }

        if (valid.n_elem > 1)
        {
          // differencing operator
          arma::sp_mat diff = arma::speye(valid.n_elem, valid.n_elem);
          for (unsigned k=0; k<_order; ++k)
          {
            arma::sp_mat d = -1.0 * arma::speye(diff.n_rows - 1, diff.n_rows);
            for(unsigned l=0; l<d.n_rows; ++l)
            {
              d.at(l,l+1) = 1.0;
            }
            diff = d * diff;
          }

          // Gaussian penalty
          // TODO: should this include normalizing constant?
          arma::vec diff_ij = diff * arma::log10(par_ij.elem(valid));
          penalty += arma::accu(-0.5 * arma::pow(diff_ij, 2) * std::pow(_penalty.at(i,j), 2));
          //penalty += arma::accu(arma::log_normpdf(diff_ij, 0, _penalty.at(i,j)));

          // gradient
          diff_ij = -1. * diff_ij * std::pow(_penalty.at(i,j), 2);
          gra_ij.elem(valid) = arma::trans(diff) * diff_ij;
          gra_ij.elem(valid) /= (par_ij.elem(valid) * log(10));
          gradient_M.tube(i,j) = gra_ij;
        }

        //OLD using all values
        //arma::vec diff_ij = diff * arma::log10(arma::vectorise(_M.tube(i,j)));
        //penalty += arma::accu(-0.5 * arma::pow(diff_ij, 2) * std::pow(_penalty.at(i,j), 2));
        //diff_ij = -1. * diff_ij * std::pow(_penalty.at(i,j), 2);
        //gradient_M.tube(i,j) = diff.t() * diff_ij; 
        //gradient_M.tube(i,j) /= (_M.tube(i,j) * log(10));
      }
    }

    return Rcpp::List::create(
        Rcpp::_["penalty"] = penalty,
        Rcpp::_["gradient"] = gradient_M
        );
  }

  Rcpp::List expected_rates (const arma::mat& _X, const arma::cube& _M, const arma::cube& _A)
  {
    /*
     *  Calculate likelihood, gradient of migration and admixture parameters
     *  given starting state
     */

    unsigned num_emissions = emission_mapping.n_cols;

    arma::mat _y = arma::zeros(num_emissions, T);
    arma::mat _B = arma::ones(num_emissions, T);

    Rcpp::List loglik = loglikelihood(_y, _B, _X, _M, _A);

    return Rcpp::List::create(
        Rcpp::_["X"] = loglik["X"], 
        Rcpp::_["y"] = loglik["y_hat"]
        );
  }

  Rcpp::List loglikelihood (const arma::mat& _y, const arma::mat& _B, arma::mat _X, const arma::cube& _M, const arma::cube& _A)
  {
    /*
     *  Calculate likelihood, gradient of migration and admixture parameters
     *  given starting state
     */

    if (_y.n_cols != T)
    {
      Rcpp::stop(prefix + "Rate statistics matrix has the wrong dimensions");
    }

    if (_B.n_cols != T)
    {
      Rcpp::stop(prefix + "Precision matrix list is the wrong dimension");
    }

    if (_M.n_slices != T)
    {
      Rcpp::stop(prefix + "Demographic parameter array has the wrong dimensions");
    }

    if (_A.n_slices != T)
    {
      Rcpp::stop(prefix + "Admixture parameter array has the wrong dimensions");
    }

    std::vector<CoalescentEpoch> epochs;
    epochs.reserve(T);

    arma::mat y_hat (emission_mapping.n_cols, T);
    arma::mat residuals (emission_mapping.n_cols, T);

    double loglik = 0.;
    for (unsigned i=0; i<T; ++i)
    {
      epochs.emplace_back(
          _X, _y.col(i), _B.col(i), _M.slice(i), _A.slice(i), t.at(i), rate_matrix_template.X,
          emission_mapping, state_mapping, initial_mapping, check_valid
      );
      y_hat.col(i) = epochs[i].y_hat;
      residuals.col(i) = epochs[i].residual;
      loglik += epochs[i].loglikelihood;
    }

    arma::mat gradient_X = arma::zeros(arma::size(_X));
    arma::cube gradient_M (arma::size(_M));
    arma::cube gradient_A (arma::size(_A));

    for (int i=T-1; i>=0; --i)
    {
      arma::cube gradient_parameters = epochs[i].reverse_differentiate(gradient_X);
      gradient_M.slice(i) = gradient_parameters.slice(0);
      gradient_A.slice(i) = gradient_parameters.slice(1);
    }

    return Rcpp::List::create(
        Rcpp::_["loglikelihood"] = loglik,
        Rcpp::_["X"] = _X, 
        Rcpp::_["y_hat"] = y_hat, 
        Rcpp::_["residuals"] = residuals,
        Rcpp::_["gradient"] = Rcpp::List::create(
          Rcpp::_["X"] = gradient_X,
          Rcpp::_["M"] = gradient_M,
          Rcpp::_["A"] = gradient_A
          )
        );
  }

  arma::mat initial_state_vectors (void)
  {
    /*
     *  Returns state vectors at time 0
     */

    arma::mat states = arma::zeros<arma::mat>(
        arma::accu(rate_matrix_template.S), 
        initial_mapping.n_cols
    );
    for (unsigned i=0; i<initial_mapping.n_cols; ++i)
    {
      arma::uword col = initial_mapping.at(0,i);
      arma::uword row = initial_mapping.at(1,i);
      states.at(row, col) = 1.0;
    }
    return states;
  }
};

//---------------------TESTS---------------------//

// [[Rcpp::export]]
Rcpp::List test_SparseMatrixExponentialMultiply (arma::sp_mat A, arma::mat B, double t, arma::mat g)
{
  SparseMatrixExponentialMultiply out (A, B, t);
  mat dB; 
  sp_mat dA = out.reverse_differentiate(dB, g);

  arma::uvec Bdim (out._Bs.size());
  for(unsigned i=0; i<out._Bs.size(); ++i) Bdim[i] = out._Bs[i].size();

  return Rcpp::List::create(
      Rcpp::_["dA"] = dA,
      Rcpp::_["dB"] = dB,
      Rcpp::_["X"] = out.result,
      Rcpp::_["Bs"] = Bdim
      );
}

// [[Rcpp::export]]
Rcpp::List test_TrioTransitionRates (arma::mat M, arma::sp_mat G)
{
  std::vector<std::string> names (M.n_rows);
  for (unsigned i=0; i<M.n_rows; ++i) names[i] = std::to_string(i);
  TrioTransitionRates rates (M);
  return Rcpp::List::create(
      Rcpp::_["matrix"] = rates.X,
      Rcpp::_["s_map_e"] = rates.states_to_emission(),
      Rcpp::_["i_map_s"] = rates.initial_to_states(),
      Rcpp::_["e_map_i"] = rates.emission_to_initial(),
      Rcpp::_["emission_states"] = rates.emission_states(names),
      Rcpp::_["initial_states"] = rates.initial_states(names),
      Rcpp::_["S"] = rates.S,
      Rcpp::_["reverse_differentiate"] = rates.reverse_differentiate(G)
  );
}

// [[Rcpp::export]]
Rcpp::List test_TrioAdmixtureProportions (arma::mat A, arma::sp_mat G)
{
  TrioAdmixtureProportions props (A, false); //disable checks for numerical differentiation
  return Rcpp::List::create(
      Rcpp::_["matrix"] = props.X,
      Rcpp::_["S"] = props.S,
      Rcpp::_["reverse_differentiate"] = props.reverse_differentiate(G)
  );
}

// [[Rcpp::export()]]
Rcpp::List test_CoalescentEpoch (arma::mat states, arma::mat M, arma::mat A, arma::mat gradient)
{
  // get maps
  TrioTransitionRates rate_matrix (M, false);
  arma::umat e_map = rate_matrix.emission_to_initial();
  arma::umat s_map = rate_matrix.states_to_emission();
  arma::umat i_map = rate_matrix.initial_to_states();
  
  // duration
  double t = 0.31;

  // get statistics
  arma::vec y = arma::ones(e_map.n_cols);
  arma::vec B = arma::ones(e_map.n_cols);

  CoalescentEpoch epoch (states, y, B, M, A, t, rate_matrix.X, e_map, s_map, i_map, false);
  arma::cube gradient_M = epoch.reverse_differentiate(gradient);

  return Rcpp::List::create(
      Rcpp::_["states_end"] = states, 
      Rcpp::_["states_start"] = epoch.states,
      Rcpp::_["gradient_states"] = gradient,
      Rcpp::_["gradient_M"] = gradient_M.slice(0),
      Rcpp::_["gradient_A"] = gradient_M.slice(1),
      Rcpp::_["y_hat"] = epoch.y_hat,
      Rcpp::_["residual"] = epoch.residual,
      Rcpp::_["gradient"] = epoch.gradient,
      Rcpp::_["loglikelihood"] = epoch.loglikelihood,
      Rcpp::_["admix_matrix"] = epoch.admix_matrix.X,
      Rcpp::_["rate_matrix"] = epoch.rate_matrix.X
      );
}

//---------------- EXPOSED CLASSES ----------------//

RCPP_MODULE(CoalescentDecoder) {
  using namespace Rcpp;
  class_<CoalescentDecoder>("CoalescentDecoder")
    .constructor<unsigned, arma::vec, bool>()
    .method("initial_states", &CoalescentDecoder::initial_states)
    .method("emission_states", &CoalescentDecoder::emission_states)
    .method("transitory_states", &CoalescentDecoder::transitory_states)
    .method("transition_rates", &CoalescentDecoder::transition_rates)
    .method("admixture_proportions", &CoalescentDecoder::admixture_proportions)
    .method("coalescence_rates", &CoalescentDecoder::coalescence_rates)
    .method("transition_operator", &CoalescentDecoder::transition_operator)
    .method("transition_operator_unsafe", &CoalescentDecoder::transition_operator_unsafe)
    .method("admixture_operator", &CoalescentDecoder::admixture_operator)
    .method("occupancy_probabilities", &CoalescentDecoder::occupancy_probabilities)
    .method("smoothness_penalty", &CoalescentDecoder::smoothness_penalty)
    .method("expected_rates", &CoalescentDecoder::expected_rates)
    .method("loglikelihood", &CoalescentDecoder::loglikelihood)
    .method("initial_state_vectors", &CoalescentDecoder::initial_state_vectors)
    ;
  class_<SparseMatrixExponentialMultiplyRescale>("SpMatExp")
    .constructor<arma::sp_mat, arma::mat, double>()
    .field("B", &SparseMatrixExponentialMultiplyRescale::_B)
    .field("Bs", &SparseMatrixExponentialMultiplyRescale::_Bs)
    .field("m_star", &SparseMatrixExponentialMultiplyRescale::_m_star)
    .field("mu", &SparseMatrixExponentialMultiplyRescale::_mu)
    .field("t", &SparseMatrixExponentialMultiplyRescale::_t)
    .field("s", &SparseMatrixExponentialMultiplyRescale::_s)
    .field("A", &SparseMatrixExponentialMultiplyRescale::_A)
    .field("result", &SparseMatrixExponentialMultiplyRescale::result)
    .method("OpNormD", &SparseMatrixExponentialMultiplyRescale::OpNormD)
    ;
  class_<OneNormEst>("OneNormEst")
    .constructor<arma::sp_mat, int, bool>()
    .field("est", &OneNormEst::_est)
    .field("w", &OneNormEst::_w)
    .field("v", &OneNormEst::_v)
    .field("iter", &OneNormEst::_iter)
    ;
}
