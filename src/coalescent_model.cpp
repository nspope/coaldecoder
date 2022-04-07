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
  const bool use_rates = false; //see note in cnstr

  const TrioAdmixtureProportions admix_matrix;
  const TrioTransitionRates rate_matrix;

  const arma::umat emission_mapping;
  const arma::umat state_mapping;
  const arma::umat initial_mapping;

  const arma::vec y; //observed statistics
  const arma::mat B; //bootstrap precision of y
  const arma::mat A; //admixture proportions
  const double t; //duration of epoch

  arma::vec y_hat;
  arma::vec residual;
  arma::vec gradient;
  arma::vec uncoalesced; //uncoalesced lineages at start
  arma::mat states; //states immediately after admixture
  double loglikelihood;

  CoalescentEpoch (
      arma::mat& _states,
      const arma::mat& _y, 
      const arma::mat& _B,
      const arma::mat& _M,
      const arma::mat& _A,
      const double& _t,
      // mappings
      const arma::umat& _emission_mapping,
      const arma::umat& _state_mapping,
      const arma::umat& _initial_mapping,
      const bool _check_valid = true
  ) : check_valid (_check_valid)
    , admix_matrix (_A, _check_valid)
    , rate_matrix (_M, _check_valid)
    , emission_mapping (_emission_mapping)
    , state_mapping (_state_mapping)
    , initial_mapping (_initial_mapping)
    , y (_y)
    , B (_B)
    , t (_t)
  {
    /* 
     *  Update states and calculate likelihood, gradient of emissions within epoch.
     */

    if (_y.n_elem != emission_mapping.n_cols)
    {
      Rcpp::stop(prefix + "Vector of rates has the wrong dimension");
    }

    if (_B.n_rows != emission_mapping.n_cols ||
        _B.n_cols != emission_mapping.n_cols )
    {
      Rcpp::stop(prefix + "Bootstrap precision matrix has the wrong dimension");
    }

    if (_states.n_rows != arma::accu(rate_matrix.S) || 
        _states.n_cols != initial_mapping.n_cols)
    {
      Rcpp::stop(prefix + "State probability vectors have the wrong dimension");
    }

    if (check_valid)
    {
      if (_t < 0.0)
      {
        Rcpp::stop(prefix + "Epoch duration cannot be negative");
      }

      if (arma::any(arma::abs(arma::sum(_states, 0) - 1.0) > _states.n_rows*arma::datum::eps))
      {
        Rcpp::stop(prefix + "State probability vector does not sum to one");
      }
    }

    unsigned num_emission = emission_mapping.n_cols;
    unsigned num_initial = initial_mapping.n_cols;

    // Admixture
    states = _states;
    _states = admix_matrix.X * _states;

    // Transition probabilities
    SparseMatrixExponentialMultiply transition (arma::trans(rate_matrix.X), _states, t);
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

    uncoalesced = arma::zeros(num_emission);
    int offset = num_emission / 2;
    for (int j=0; j<num_emission; ++j)
    {
      arma::uword col = emission_mapping.at(0,j);
      arma::uword lin = emission_mapping.at(1,j);

      // 2-lineage states can transition to 1-lineage states
      uncoalesced[j] = lin == 2 ? transitory_start[col] :
        transitory_start[col] + coalesced_start[j - offset]; 

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
    if (use_rates)
    {
      y_hat /= t;
    }

    // Calculate loglikelihood
    residual = y - y_hat;
    gradient = B * residual;
    loglikelihood = -0.5 * arma::dot(residual, gradient);
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

    // Gradient wrt state vectors
    arma::mat d_states = arma::zeros(arma::size(states));
    arma::vec d_uncoalesced = -y_hat % gradient/uncoalesced;
    arma::vec d_coalesced_start = -gradient/uncoalesced;
    arma::vec d_coalesced_end = gradient/uncoalesced;
    if (use_rates)
    {
      // TODO check this
      d_coalesced_start /= t;
      d_coalesced_end /= t;
    }
    arma::vec d_transitory_start = arma::zeros(num_initial);
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
        d_coalesced_start[j - offset] += d_uncoalesced[j];
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
    SparseMatrixExponentialMultiply transition (
      arma::trans(rate_matrix.X), admix_matrix.X * states, t);
    arma::mat tmp = arma::zeros(arma::size(_states));
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
    //arma::sp_mat d_admix_matrix (_states * arma::trans(states));
    _states = arma::trans(admix_matrix.X) * _states;
    _states = _states + d_states;

    // Gradient wrt demographic parameter, admixture matrices
    arma::cube d_parameters (rate_matrix.P, rate_matrix.P, 2);
    d_parameters.slice(0) = rate_matrix.reverse_differentiate(arma::trans(d_rate_matrix));
    d_parameters.slice(1) = admix_matrix.reverse_differentiate(d_admix_matrix);

    return d_parameters;
  }

  //arma::cube forward_differentiate (arma::mat& _states)
  //{
  //  // want to calculate hessian if possible
  //  // input: dll / dpar_i
  //  // want to get: (dll / dpar_i) / dpar_j
  //  // the tricky part is the matrix exponential
  //  // in the forward pass:
  //  //   parameters -> rate matrix -> matrix exponential multiply
  //  //   drate_matrix/dparameter * dmatrix exponential multiply/drate_matrix
  //  // so we need to be able to multiply by the *adjoint* of the jacobian for the expmat
  //  // we'd also need to differentiate the jacobian of the expmet wrt other rate parameters
  //  // this seems hard, let's pass for now
  //  // it might be possible to use an autodiff library
  //}
};

struct CoalescentDecoder
{
  const std::string prefix = "[CoalescentDecoder] ";
  const bool check_valid = true;
  const bool use_marginal_statistics = true; //see notes in method coalescence_rate
  const bool use_rates = false; //divide coalescence probabilities by time

  const unsigned P; //number of populations
  const unsigned T; //number of epochs
  const TrioTransitionRates rate_matrix_template;
  const arma::umat emission_mapping;
  const arma::umat state_mapping;
  const arma::umat initial_mapping;

  arma::mat y; // input statistics
  arma::cube B; // bootstrap precision
  arma::vec t; // epoch duration

  // constructor with no input data
  CoalescentDecoder (
      const unsigned _P, 
      const arma::vec& _t, // epoch duration
      const bool _use_rates = true,
      const bool _check_valid = true
  ) : check_valid (_check_valid)
    , use_rates (_use_rates)
    , P (_P)
    , T (_t.n_elem)
    , rate_matrix_template (arma::ones(_P, _P))
    , emission_mapping (rate_matrix_template.emission_to_initial())
    , state_mapping (rate_matrix_template.states_to_emission())
    , initial_mapping (rate_matrix_template.initial_to_states())
  {
    t = _t;
    y = arma::mat(emission_mapping.n_cols, T, arma::fill::zeros);
    B = arma::cube(y.n_rows, y.n_rows, T, arma::fill::zeros);
  }

  // constructor with input data
  CoalescentDecoder (
      const unsigned _P, 
      const arma::cube& _y, // number lineages coalescing in interval
      const arma::mat& _n, // number lineages 
      const arma::vec& _t, // epoch duration
      const bool _use_rates = true,
      const bool _check_valid = true
  ) : check_valid (_check_valid)
    , use_rates (_use_rates)
    , P (_P)
    , T (_t.n_elem)
    , rate_matrix_template (arma::ones(_P, _P))
    , emission_mapping (rate_matrix_template.emission_to_initial())
    , state_mapping (rate_matrix_template.states_to_emission())
    , initial_mapping (rate_matrix_template.initial_to_states())
  {
    if (_y.n_rows != emission_mapping.n_cols || _y.n_cols != T)
    {
      Rcpp::stop(prefix + " dimension of 'y' does not match");
    }
    if (_n.n_rows != emission_mapping.n_cols || _n.n_cols != _y.n_slices)
    {
      Rcpp::stop(prefix + " dimension of 'n' does not match");
    }

    // epoch durations
    t = _t;

    // convert raw statistics to rates
    arma::cube x (arma::size(_y));
    for (unsigned i=0; i<_y.n_slices; ++i)
    {
      x.slice(i) = coalescence_rates(_y.slice(i), _n.col(i), _t, use_rates, use_marginal_statistics);
    }

    // observed statistics
    y = x.slice(0);

    // bootstrap precision matrices
    // (use identity matrix if no bootstrap reps)
    //
    // TODO: re-evaluate strategy here. If inputs are small (they are)
    // this should be internally rescaled. If there are a lot of pops, we may want
    // to use a shrinkage approach, which would require implementing something separately.
    // As a workaround for now, maybe add a separate constructor that allows B to be specified.
    B = arma::cube(_y.n_rows, _y.n_rows, T, arma::fill::zeros);
    unsigned num_boot = _y.n_slices - 1;
    unsigned offset = emission_mapping.n_cols / 2;
    for (unsigned i=0; i<T; ++i)
    {
      if (num_boot > 0)
      {
        arma::mat b = arma::zeros(num_boot, x.n_rows);
        for (unsigned j = 0; j < num_boot; ++j)
        {
          b.row(j) = arma::trans(x.slice(j+1).col(i));
        }
        //<DEBUG> "disable" trio rates
        auto pair_idx = arma::span(0, offset-1);//DEBUG
        b = b.cols(pair_idx);//DEBUG
        B.slice(i).submat(pair_idx, pair_idx) = arma::inv(arma::cov(b));//DEBUG
        //<\DEBUG>
        //B.slice(i) = arma::inv(arma::cov(b));
      } else {
        //<DEBUG> "disable" trio rates
        auto pair_idx = arma::span(0, offset-1);//DEBUG
        B.slice(i).submat(pair_idx, pair_idx).eye();
        //<\DEBUG>
        //B.slice(i).eye();
      }
    }
  }

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

  arma::mat observed_rates (void) const
  {
    return y;
  }

  arma::cube precision_matrices (void) const
  {
    return B;
  }

  std::vector<std::string> initial_states (void) const
  {
    return rate_matrix_template.initial_states();
  }

  std::vector<std::string> emission_states (void) const
  {
    return rate_matrix_template.emission_states();
  }

  arma::sp_mat transition_rates (const arma::mat& _M)
  {
    /*
     *  Return trio transition rate matrix given demographic parameters
     */

    TrioTransitionRates rates (_M, check_valid);
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

  arma::mat transition_operator (const mat& _X, const arma::mat& _M, const double& _t)
  {
    /*
     *  Map state probability vectors across duration, given demographic
     *  parameters
     */
    TrioTransitionRates rates (_M, check_valid);
    SparseMatrixExponentialMultiply expm (arma::trans(rates.X), _X, _t);
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
      arma::vec Asum = arma::sum(_A, 1);
      if (arma::any(arma::abs(Asum - 1.0) > _A.n_cols*arma::datum::eps))
      {
        Rcpp::stop(prefix + "Admixture proportions do not sum to one");
      }
    }

    arma::cube out = arma::zeros<arma::cube>(P, P, T + 1);
    out.slice(0).eye();
    for (unsigned i=0; i<T; ++i)
    {
      arma::mat rate_matrix = _M.slice(i);
      for (unsigned j=0; j<P; ++j)
      {
        rate_matrix.at(i,i) = 0.;
        rate_matrix.at(i,i) = -arma::accu(rate_matrix.row(i));
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

    // construct difference operator
    arma::sp_mat diff = arma::speye(T, T);
    for (unsigned i=0; i<_order; ++i)
    {
      arma::sp_mat d = -1.0 * arma::speye(diff.n_rows - 1, diff.n_rows);
      for(unsigned j=0; j<d.n_rows; ++j)
      {
        d.at(j,j+1) = 1.0;
      }
      diff = d * diff;
    }

    // prior, gradient of prior wrt parameters
    double penalty = 0;
    for (unsigned i=0; i<P; ++i)
    {
      for (unsigned j=0; j<P; ++j)
      {
        arma::vec diff_ij = diff * arma::log10(arma::vectorise(_M.tube(i,j)));
        penalty += arma::accu(-0.5 * arma::pow(diff_ij, 2) * std::pow(_penalty.at(i,j), 2));
        diff_ij = -1. * diff_ij * std::pow(_penalty.at(i,j), 2);
        gradient_M.tube(i,j) = diff.t() * diff_ij; 
        gradient_M.tube(i,j) /= (_M.tube(i,j) * log(10));
      }
    }

    return Rcpp::List::create(
        Rcpp::_["penalty"] = penalty,
        Rcpp::_["gradient"] = gradient_M
        );
  }

  Rcpp::List loglikelihood (arma::mat _X, const arma::cube& _M, const arma::cube& _A)
  {
    /*
     *  Calculate likelihood, gradient of migration and admixture parameters
     *  given starting state
     */

    if (_M.n_slices != T)
    {
      Rcpp::stop(prefix + "Demographic parameter array has wrong dimensions");
    }

    if (_A.n_slices != T)
    {
      Rcpp::stop(prefix + "Admixture parameter array has wrong dimensions");
    }

    std::vector<CoalescentEpoch> epochs;
    epochs.reserve(T);

    arma::mat y_hat (emission_mapping.n_cols, T);
    arma::mat residuals (emission_mapping.n_cols, T);

    double loglik = 0.;
    for (unsigned i=0; i<T; ++i)
    {
      epochs.emplace_back(
          _X, y.col(i), B.slice(i), _M.slice(i), _A.slice(i), t.at(i), 
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
  TrioTransitionRates rates (M);
  return Rcpp::List::create(
      Rcpp::_["matrix"] = rates.X,
      Rcpp::_["s_map_e"] = rates.states_to_emission(),
      Rcpp::_["i_map_s"] = rates.initial_to_states(),
      Rcpp::_["e_map_i"] = rates.emission_to_initial(),
      Rcpp::_["emission_states"] = rates.emission_states(),
      Rcpp::_["initial_states"] = rates.initial_states(),
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
  arma::mat P = arma::eye(e_map.n_cols, e_map.n_cols);

  CoalescentEpoch epoch (states, y, P, M, A, t, e_map, s_map, i_map, false);
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
    .constructor<unsigned, arma::cube, arma::mat, arma::vec, bool, bool>()
    .constructor<unsigned, arma::vec, bool, bool>()
    .method("observed_rates", &CoalescentDecoder::observed_rates)
    .method("precision_matrices", &CoalescentDecoder::precision_matrices)
    .method("initial_states", &CoalescentDecoder::initial_states)
    .method("emission_states", &CoalescentDecoder::emission_states)
    .method("transition_rates", &CoalescentDecoder::transition_rates)
    .method("admixture_proportions", &CoalescentDecoder::admixture_proportions)
    .method("coalescence_rates", &CoalescentDecoder::coalescence_rates)
    .method("transition_operator", &CoalescentDecoder::transition_operator)
    .method("admixture_operator", &CoalescentDecoder::admixture_operator)
    .method("occupancy_probabilities", &CoalescentDecoder::occupancy_probabilities)
    .method("smoothness_penalty", &CoalescentDecoder::smoothness_penalty)
    .method("loglikelihood", &CoalescentDecoder::loglikelihood)
    .method("initial_state_vectors", &CoalescentDecoder::initial_state_vectors)
    ;
}
