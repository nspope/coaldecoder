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

      if (arma::any(arma::sum(_states, 0) != 1.0))
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
    arma::sp_mat d_rate_matrix = 
      transition.reverse_differentiate(d_states, _states); 
    _states = d_states;

    // Gradient wrt pre-admixture state vectors, trio admixture matrix
    arma::sp_mat d_admix_matrix = admix_matrix.X;
    for (arma::sp_mat::iterator it = d_admix_matrix.begin(); 
         it != d_admix_matrix.end(); ++it)
    {
      double val = (*it) + arma::dot(_states.row(it.row()), states.row(it.col()));
      if (val == 0.0)
      {
        Rcpp::stop(prefix + " Can't void nonzero element");
      }
      (*it) = val;
    }
    d_admix_matrix -= admix_matrix.X;
    _states = arma::trans(admix_matrix.X) * _states;

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
  const bool use_marginal_statistics = true;

  const unsigned P; //number of populations
  const unsigned T; //number of epochs
  const TrioTransitionRates rate_matrix_template;
  const arma::umat emission_mapping;
  const arma::umat state_mapping;
  const arma::umat initial_mapping;

  arma::mat y; // input statistics
  arma::cube B; // bootstrap precision
  arma::vec t; // epoch duration

  CoalescentDecoder (
      const unsigned _P, 
      const arma::mat& _y, 
      const arma::cube& _y_boot, 
      const arma::vec& _t,
      const bool _check_valid
  ) : check_valid (_check_valid)
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
    if (_y_boot.n_rows != emission_mapping.n_cols || _y_boot.n_cols != T)
    {
      Rcpp::stop(prefix + " dimension of 'y_boot' does not match");
    }

    // process data, estimate bootstrap precision matrices
    y = _y;
    B = arma::cube(_y.n_rows, _y.n_rows, T);
    t = _t;
    for (unsigned i=0; i<T; ++i)
    {
      arma::mat x = _y_boot.col(i);
      B.slice(i) = arma::inv(arma::cov(arma::trans(x)));
    }
    if (!use_marginal_statistics)
    {
      //TODO
      // Input rates are marginal:
      //   [p(t > end) - p(t > start)] / [1.0 - p(t < start)]
      // We want:
      //   [p(t1 > end & t2 > end) - p(t1 > start & t2 > end)] / [1.0 - p(t1 < start)]
      // To calculate these we need to modify numerator and denominator
    }
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
      if (arma::any(Asum != 1.0))
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

  Rcpp::List smoothness_penalty (const arma::cube& _M, const arma::mat& penalty, const unsigned order = 1)
  {
    /*
     *  Calculates a penalty on squared differences of a given order, for each
     *  demographic parameter. The differencing is done with respect to time
     *  and on a log10 scale.  Also calculates the gradient wrt the penalty.
     */

    if (order == 0)
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
    if (penalty.n_rows != P || penalty.n_cols != P)
    {
      Rcpp::stop(prefix + "Penalty matrix has wrong dimensions");
    }
    if (arma::any(arma::vectorise(penalty) < 0.0))
    {
      Rcpp::stop(prefix + "Penalty matrix has negative values");
    }

    arma::cube gradient_M (arma::size(_M), arma::fill::zeros);

    // construct difference operator
    arma::sp_mat diff = arma::speye(T, T);
    for (unsigned i=0; i<order; ++i)
    {
      arma::sp_mat d = -1.0 * arma::speye(diff.n_rows - 1, diff.n_rows);
      for(unsigned j=0; j<d.n_rows; ++j)
      {
        d.at(j,j+1) = 1.0;
      }
      diff = d * diff;
    }

    // prior, gradient of prior wrt parameters
    double logprior = 0;
    for (unsigned i=0; i<P; ++i)
    {
      for (unsigned j=0; j<P; ++j)
      {
        arma::vec diff_ij = diff * arma::log10(arma::vectorise(_M.tube(i,j)));
        logprior += arma::accu(-0.5 * arma::pow(diff_ij, 2) * std::pow(penalty.at(i,j), 2));
        diff_ij = -1. * diff_ij * std::pow(penalty.at(i,j), 2);
        gradient_M.tube(i,j) = diff.t() * diff_ij; 
        gradient_M.tube(i,j) /= (_M.tube(i,j) * log(10));
      }
    }

    return Rcpp::List::create(
        Rcpp::_["logprior"] = logprior,
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

  arma::mat initial_states (void)
  {
    /*
     *  Returns state vectors at time 0
     */

    arma::mat states = arma::zeros<arma::mat>(
        arma::accu(rate_matrix_template.S), 
        arma::accu(rate_matrix_template.S)
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
  double ll = epoch.loglikelihood;
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

