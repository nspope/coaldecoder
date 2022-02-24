#include <RcppArmadillo.h> 
#include <vector>
#include <nloptrAPI.h>
#include "expm.h"
#include "expm_revdiff.h"

// [[Rcpp::export]]
Rcpp::List test_matrix_exponential_multiply (arma::mat A, arma::mat B, double t, arma::mat g)
{
  matrix_exponential_multiply out (A, B, t);
  mat dA; mat dB; 
  out.reverse_differentiate(dA, dB, g);
  return Rcpp::List::create(
      Rcpp::_["dA"] = dA,
      Rcpp::_["dB"] = dB,
      Rcpp::_["X"] = out.result
      );
}

struct transition_rate_matrix
{
  public:
  unsigned P, I;
  arma::umat t_states, c_states, e_states, s_states, c_map_e, s_mask_e;
  arma::uvec s_map_t;
  unsigned num_trans, num_coal, num_start, num_emiss;

  const bool approx; // use sparser approximation where state transitions can involve only one migration

  transition_rate_matrix (const unsigned P, const unsigned I, const bool approx = true) : P(P), I(I), approx(approx)
  {
    // enumerate transitory states (e.g. haplotype configurations across populations)
    arma::uvec pop_index = arma::regspace<arma::uvec>(0, P-1);
    t_states = pop_index;
    for (unsigned i=1; i<I; ++i)
    {
      t_states = arma::join_horiz(arma::repmat(t_states, pop_index.n_elem, 1), arma::repelem(pop_index, t_states.n_rows, 1));
    }

    // find nonredundant starting configurations
    s_map_t = arma::regspace<arma::uvec>(0, t_states.n_rows-1);
    std::sort(s_map_t.begin(), s_map_t.end(), [&](unsigned i, unsigned j){
              arma::urowvec a = arma::sort(t_states.row(i)),
                            b = arma::sort(t_states.row(j));
              for(unsigned k=0; k<I-1; ++k) 
              {
                if(a.at(k) != b.at(k))
                {
                   return a.at(k) < b.at(k);
                }
              }
              return a.at(I-1) < b.at(I-1);
    });
    auto end = std::unique(s_map_t.begin(), s_map_t.end(), [&](unsigned i, unsigned j){
                arma::irowvec a = arma::conv_to<arma::irowvec>::from(arma::sort(t_states.row(i))),
                              b = arma::conv_to<arma::irowvec>::from(arma::sort(t_states.row(j)));
                if (arma::any(a - b))
                  return false;
                return true;
    });
    s_map_t = s_map_t.head(end - s_map_t.begin());
    s_states = arma::umat(s_map_t.n_elem, P, arma::fill::zeros);
    for (unsigned i=0; i<s_map_t.n_elem; ++i)
    {
      for (unsigned j=0; j<I; ++j)
      {
        s_states.at(i,t_states.at(s_map_t.at(i),j)) += 1;
      }
    }

    // enumerate possible pairwise coalescence events
    c_states = arma::zeros<arma::umat>(0, 2);
    for (unsigned i=0; i<I; ++i)
    {
      for (unsigned j=i; j<I; ++j)
      {
        if (i != j)
        {
          c_states = arma::join_vert(c_states, arma::urowvec({i,j}));
        }
      }
    }

    // find nonredundant coalescent events for a given starting state
    // e.g. conditional on starting state, divide coal events between indexed lineages into within/cross population
    e_states = arma::zeros<arma::umat>(0, 2);
    for (unsigned i=0; i<P; ++i)
    {
      for (unsigned j=i; j<P; ++j)
      {
        e_states = arma::join_vert(e_states, arma::urowvec({i,j}));
      }
    }
    c_map_e = arma::zeros<arma::umat>(c_states.n_rows, s_map_t.n_elem);
    s_mask_e = arma::ones<arma::umat>(e_states.n_rows, s_map_t.n_elem);
    for(unsigned i=0; i<s_map_t.n_elem; ++i)
    {
      for(unsigned j=0; j<c_states.n_rows; ++j)
      {
        unsigned p = t_states.at(s_map_t.at(i),c_states.at(j,0)), 
                 q = t_states.at(s_map_t.at(i),c_states.at(j,1));
        arma::uvec idx = arma::find(e_states.col(0) == std::min(p,q) && e_states.col(1) == std::max(p,q), 1);
        if (idx.n_elem)
        {
          c_map_e.at(j,i) = idx.at(0);
          s_mask_e.at(c_map_e.at(j,i),i) = 0;
        }
        else
        {
          Rcpp::warning("coalescent event does not map to emission");
        }
      }
    }

    num_trans = t_states.n_rows;
    num_coal  = c_states.n_rows;
    num_start = s_states.n_rows;
    num_emiss = e_states.n_rows;
  }

  std::string symbolic_transition_rates (const arma::uvec& remap)
  {
    if (!(remap.n_elem == P && arma::max(remap) < P)) Rcpp::stop("invalid remap");

    std::string out = "";
    std::string outstring;

    // transition rates to transient states
    for (unsigned i=0; i<num_trans; ++i)
      for (unsigned j=0; j<num_trans; ++j)
      {
        if (i != j)
        {
          unsigned moves = 0;
          outstring  = "M[" + std::to_string(i+1) + "," + std::to_string(j+1) + "] = 1";
          for (unsigned k=0; k<I; ++k)
          {
            unsigned p1 = remap.at(t_states.at(i,k)), // origin, backwards in time; destination, forwards in time
                     p2 = remap.at(t_states.at(j,k)); // destination, backwards in time; origin, forwards in time
            if (p1 != p2)
            {
              moves += 1;
              outstring += "*m[" + std::to_string(p1+1) + "," + std::to_string(p2+1) + "]";
            }
          }
          if (moves == 0) 
            outstring = "M[" + std::to_string(i+1) + "," + std::to_string(j+1) + "] = 0";
          outstring += ";\n";
          out += outstring;
        }
      }

    // transition rates to coalesced states
    for (unsigned i=0; i<num_trans; ++i)
    {
      for (unsigned j=0; j<num_coal; ++j)
      {
        const unsigned p1 = remap.at(t_states.at(i,c_states.at(j,0))),
                       p2 = remap.at(t_states.at(i,c_states.at(j,1)));
        outstring = "M[" + std::to_string(i+1) + "," + std::to_string(num_trans+j+1) + "]";
        outstring += p1==p2 ? 
          " = 1/m[" + std::to_string(p1+1) + "," + std::to_string(p1+1) + "];\n" :
          " = 0;\n";
        out += outstring;
      }
    }

    // total exit rates
    for (unsigned i=0; i<num_trans+num_coal; ++i)
    {
      outstring = "M[" + std::to_string(i+1) + "," + std::to_string(i+1) + "] = ";
      for (unsigned j=0; j<num_trans+num_coal; ++j)
        if (i != j)
          outstring += "-M[" + std::to_string(i+1) + "," + std::to_string(j+1) + "]";
      outstring += ";\n";
      out += outstring;
    }

   return out;
  }

  arma::mat transition_rates (const arma::uvec& remap, arma::mat migr_mat)
  {
    // Transition rate matrix for migration process up until first coalescent event
    // with populations remapped to a population vector of lower order. This results in a transition
    // matrix with redundant states. 

    arma::vec ne = migr_mat.diag();
    migr_mat.diag().ones();

    if (!(remap.n_elem == P && arma::max(remap) < P)) Rcpp::stop("invalid remap");
    if (!(migr_mat.n_rows == P && migr_mat.n_cols == P && arma::all(arma::vectorise(migr_mat) >= 0.))) Rcpp::stop("invalid migr_mat");
    if (!(arma::all(ne >= 0.))) Rcpp::stop("invalid ne");

    arma::mat rates (num_trans+num_coal, num_trans+num_coal, arma::fill::zeros);

    // transition rates to transient states
    for (unsigned i=0; i<num_trans; ++i)
    {
      for (unsigned j=0; j<num_trans; ++j)
      {
        unsigned moves = 0;
        rates.at(i,j) = 1.;
        for (unsigned k=0; k<I; ++k)
        {
          unsigned p1 = remap.at(t_states.at(i,k)), // origin, backwards in time; destination, forwards in time
                   p2 = remap.at(t_states.at(j,k)); // destination, backwards in time; origin, forwards in time
          if (p1 != p2)
          {
            moves += 1;
            rates.at(i,j) *= migr_mat.at(p1,p2);
          }
        }
        if (moves == 0) rates.at(i,j) = 0.;
        if (approx && moves > 1) rates.at(i,j) = 0.;
      }
    }

    // transition rates to coalesced states
    for (unsigned i=0; i<num_trans; ++i)
    {
      for (unsigned j=0; j<num_coal; ++j)
      {
        const unsigned p1 = remap.at(t_states.at(i,c_states.at(j,0))),
                       p2 = remap.at(t_states.at(i,c_states.at(j,1)));
        rates.at(i,num_trans+j) = p1==p2 ? 1./ne[p1] : 0.;
      }
    }

    // total exit rates
    for (unsigned i=0; i<rates.n_rows; ++i)
    {
      rates.at(i,i) = 0.;
      for (unsigned j=0; j<rates.n_cols; ++j)
      {
        if (i != j)
        {
          rates.at(i,i) -= rates.at(i,j);
        }
      }
    }

   return rates;
  }

  arma::mat reverse_differentiate (arma::mat d_rates, const arma::mat& rates, const arma::uvec& remap, arma::mat migr_mat)
  {
    //would be nicer to have this return a lambda with redundant inputs captured

    // Transition rate matrix for migration process up until first coalescent event
    // with populations remapped to a population vector of lower order. This results in a transition
    // matrix with redundant states. 

    arma::vec ne = migr_mat.diag();
    migr_mat.diag().ones();

    if (!(remap.n_elem == P && arma::max(remap) < P)) Rcpp::stop("invalid remap");
    if (!(migr_mat.n_rows == P && migr_mat.n_cols == P && arma::all(arma::vectorise(migr_mat) >= 0.))) Rcpp::stop("invalid migr_mat");
    if (!(arma::all(ne >= 0.))) Rcpp::stop("invalid ne");
    if (!(rates.n_cols == d_rates.n_cols && rates.n_rows == d_rates.n_rows && d_rates.n_cols == d_rates.n_rows)) Rcpp::stop("invalid rates");
    if (!(rates.n_cols == num_trans + num_coal)) Rcpp::stop("invalid rates");

    // RD: total exit rates
    for (int i=int(rates.n_rows)-1; i>=0; --i)
    {
      for (int j=int(rates.n_cols)-1; j>=0; --j)
      {
        if (i != j)
        {
          d_rates.at(i,j) -= d_rates.at(i,i);
        }
      }
      d_rates.at(i,i) = 0.;
    }

    // RD: transition rates to coalesced states
    arma::vec d_ne = arma::zeros<arma::vec>(arma::size(ne));
    for (int i=int(num_trans)-1; i>=0; --i)
    {
      for (int j=int(num_coal)-1; j>=0; --j)
      {
        const unsigned p1 = remap.at(t_states.at(i,c_states.at(j,0))),
                       p2 = remap.at(t_states.at(i,c_states.at(j,1)));
        if (p1==p2) d_ne[p1] += -d_rates.at(i,num_trans+j)/std::pow(ne[p1], 2);
      }
    }

    // RD: transition rates to transient states
    arma::mat d_migr_mat = arma::zeros(arma::size(migr_mat));
    arma::mat tmp_d (P, P);
    for (int i=int(num_trans)-1; i>=0; --i)
    {
      for (int j=int(num_trans)-1; j>=0; --j)
      {
        unsigned moves = 0;
        tmp_d.zeros();
        for (unsigned k=0; k<I; ++k)
        {
          unsigned p1 = remap.at(t_states.at(i,k)), // origin, backwards in time; destination, forwards in time
                   p2 = remap.at(t_states.at(j,k)); // destination, backwards in time; origin, forwards in time
          if (p1 != p2)
          {
            moves += 1;
            tmp_d.at(p1,p2) += d_rates.at(i,j) * rates.at(i,j)/migr_mat.at(p1,p2);
          }
        }
        //if (moves == 0) tmp_d.zeros(); //unnecessary
        if (approx && moves > 1) tmp_d.zeros();
        d_migr_mat += tmp_d;
        //for (unsigned k=0; k<I; ++k)
        //{
        //  unsigned p1 = remap.at(t_states.at(i,k)), // origin, backwards in time; destination, forwards in time
        //           p2 = remap.at(t_states.at(j,k)); // destination, backwards in time; origin, forwards in time
        //  if (p1 != p2)
        //  {
        //    d_migr_mat.at(p1,p2) += d_rates.at(i,j) * rates.at(i,j)/migr_mat.at(p1,p2);
        //  }
        //}
      }
    }

    d_migr_mat.diag() = d_ne;

    return d_migr_mat;
  }
};

//struct transition_rate_matrix_full
//{
//  // WIP
//  public:
//  unsigned P;
//  arma::umat t_states, c_states, e_states, s_states, c_map_e, s_mask_e;
//  arma::uvec s_map_t;
//  unsigned num_trans, num_coal, num_start, num_emiss;
//
//  const bool approx; // use sparser approximation where state transitions can involve only one migration
//
//  transition_rate_matrix_full (const unsigned P, const unsigned I, const bool approx = true) : P(P), approx(approx)
//  {}
//    // would have to find all states that lead to a given state
//    // so if state is ivec
//    // population index or negative
//    // negative means coalesced
//    //
//    // i can see a way to do this easily with 3
//    //
//    // our state vector is four integers
//    // first three are tips
//    // fourth is internal node
//    // 0 means the node hasnt been created yet
//    // negative means the node coalesced with whatever this number is
//    // eg
//    // -2 -1 3 1
//    // means that node 1 coalesced with 2, node 2 coalesced with 1, node 3 is in population 3, node 4 is in population 1
//    //
//    // we need to find all states that lead to the current state
//    // so: 
//    // if last node is negative, find partner of last node. freeze other two. loop over populations. parents are <frozen> <pop> <pop>
//    // if last node is positive, loop over indices. 
//    //   if negative, find partner, freeze others. loop over pops. parents are <pop> <pop> <frozen>
//    //   if positive, freeze others. loop over pops. parents are <frozen> <pop>
//    // if last node is zero, loop over indices
//    //   if negative, something went horribly wrong
//    //   if positive, freeze others. loop over pops. parents are <frozen> <pop>
//
//  arma::vec loop_over_parent_nodes (const arma::ivec& state, const arma::vec& input_vector)
//  {
//    if (input_vector.n_elem != _dim_) Rcpp::stop("input has incorrect dimension");
//    arma::vec output_vector = arma::zeros(_dim_);
//
//    output_vector[state] = _total_exit_rate_.at(state);
//
//    if (state.at(3) < 0)
//    {
//      // we're in an absorbing state
//      arma::ivec prior_state = state;
//      int mate = abs(state.at(3)) - 1;
//      for (unsigned p=0; p<P; ++p)
//      {
//        prior_state[mate] = prior_state[3] = p;
//        output_vector[_state_index_[prior_state]] += 
//          _coalescence_rate_.at(prior_state[3]) * 
//          input_vector[_state_index_[prior_state]];
//      }
//    } 
//    if (state.at(3) > 0)
//    {
//      // we're after the first coalescent event
//      for (unsigned i=0; i<4; ++i)
//      {
//        arma::ivec prior_state = state;
//        if (state.at(i) > 0)
//        {
//          for (unsigned p=0; p<P; ++p)
//          {
//            prior_state[i] = p;
//            output_vector[_state_index_[prior_state]] += 
//              _backwards_migration_from_to_.at(prior_state[i], state[i]) * 
//              input_vector[_state_index_[prior_state]];
//          }
//        }
//        if (state.at(i) < 0)
//        {
//          int mate = abs(state.at(3)) - 1;
//          for (unsigned p=0; p<P; ++p)
//          {
//            prior_state[mate] = prior_state[i] = p;
//            output_vector[_state_index_[prior_state]] += 
//              _coalescence_rate_.at(prior_state[i]) * 0.5 *
//              input_vector[_state_index_[prior_state]];
//          }
//        }
//      }
//    }
//    if (state.at(3) == 0)
//    {
//      // we're before the first coalescent event
//      for (unsigned i=0; i<3; ++i)
//      {
//        arma::ivec prior_state = state;
//        for (unsigned p=0; p<P; ++p)
//        {
//          if (i != p)
//          {
//            prior_state[i] = p;
//            output_vector[_state_index_[prior_state]] += 
//              _backwards_migration_from_to_.at(prior_state[i], state[i]) * 
//              input_vector[_state_index_[prior_state]];
//          }
//        }
//      }
//    }
//
//    return output_vector;
//  }
//};

struct penalized_migration_function
{
  transition_rate_matrix rate_matrix;

  unsigned P;
  arma::uvec remap;
  arma::mat migr_mat, migr_fun;
  double duration, penalty = 0., loglikelihood = 0.;
  bool normalizing_constant = false;

  matrix_exponential_multiply transition;

  penalized_migration_function (arma::mat& states, arma::uvec remap_in, arma::mat migr_mat_in, double duration_in, double penalty_in)
    : rate_matrix (migr_mat_in.n_rows, 1)
    , P (migr_mat_in.n_rows) 
    , remap (remap_in)
    , migr_mat (migr_mat_in)
    , duration (duration_in)
    , penalty (penalty_in)
    , loglikelihood (0.)
    , transition (arma::trans(rate_matrix.transition_rates(remap, migr_mat)), states, duration)
  {
    //if(!(states.n_rows == P && states.n_cols == P && arma::all(arma::vectorise(states) >= 0.))) Rcpp::stop("invalid states");
    if(!(states.n_rows == P && states.n_cols == P)) Rcpp::stop("invalid states");

    states = transition.result;

    migr_fun = states; //member

    // mixture-of-Dirichlet sparsity penalty:
    //   \sum_i \theta_i^penalty = \sum_i \theta_i^{penalty + 1 - 1} \prod_{j \neq i} \theta_j^{1 - 1}
    // which is the sum of kernels of a Dirichlet.
    // the normalizing constant (shared by components) is:
    //   \gamma(penalty + n - 1) / \gamma(penalty) \times 1 / \gamma(1)^(n-1) = \gamma(penalty + n - 1) / \gamma(penalty)
    // the mixture weights are 1/n
    double lp = 0., constant = 0.;
    if (penalty > 0.)
    {
      for (unsigned i=0; i<P; ++i)
      {
        // safe log sum exp
        arma::vec log_states = penalty * arma::log(migr_fun.col(i));
        double log_max = log_states.max();
        double log_sum = log_max + log(arma::accu(arma::exp(log_states - log_max)));
        lp += log_sum;
        if (normalizing_constant) constant += 1. * ::Rf_lgammafn(P - 1 + penalty) - 1. * ::Rf_lgammafn(penalty) - log(double(P));
      }
    }
    lp += constant;

    loglikelihood = lp; //member
  }

  arma::mat reverse_differentiate (arma::mat& d_states, const arma::mat& d_penalty) 
  {
    // returns gradient of rate matrix with regard to penalty function
    // updates gradient of state vectors with regard to penalty function

    if(!(d_states.n_rows == P && d_states.n_cols == P)) Rcpp::stop("invalid d_states");

    // d \log(\sum \theta^penalty) / d \theta_i = penalty \theta^{penalty - 1} / \sum \theta^penalty
    if (penalty > 0.)
    {
      for (unsigned i=0; i<P; ++i)
      {
        arma::vec log_states = penalty * arma::log(migr_fun.col(i));
        double log_max = log_states.max();
        double log_sum = log_max + log(arma::accu(arma::exp(log_states - log_max)));
        arma::vec tmp = arma::exp(log(penalty) + (penalty - 1.) * arma::log(migr_fun.col(i)) - log_sum);
        d_states.col(i) += tmp;
      }
    }

    d_states += d_penalty;//derivative of dirichlet penalty wrt states at current timepoint

    // current timepoint --> previous timepoint
    arma::mat d_states_updated, d_rates;
    transition.reverse_differentiate(d_rates, d_states_updated, d_states);
    d_states = d_states_updated;

    return rate_matrix.reverse_differentiate(arma::trans(d_rates), arma::trans(transition()), remap, migr_mat);
  }
};

struct coalescent_epoch
{
  transition_rate_matrix* rate_matrix;

  unsigned P, I, num_trans, num_coal, num_start, num_emiss;
  arma::umat c_map_e, s_mask_e, y;
  arma::uvec z, remap;
  arma::mat migr_mat, tran_probs, coal_probs, emiss_probs;
  double duration, loglikelihood;
  bool normalizing_constant = false;

  matrix_exponential_multiply transition;

  coalescent_epoch (arma::mat& states, arma::uvec& n, arma::umat& y_in, arma::uvec remap_in, arma::mat migr_mat_in, double duration_in, transition_rate_matrix* rate_matrix, bool simulate = false)
    : rate_matrix (rate_matrix)
    , P (rate_matrix->P) 
    , I (rate_matrix->I)
    , num_trans (rate_matrix->num_trans)
    , num_coal (rate_matrix->num_coal)
    , num_start (rate_matrix->num_start)
    , num_emiss (rate_matrix->num_emiss)
    , c_map_e (rate_matrix->c_map_e)
    , s_mask_e (rate_matrix->s_mask_e)
    , y (y_in)
    , remap (remap_in)
    , migr_mat (migr_mat_in)
    , duration (duration_in)
    , loglikelihood (0.)
    , transition (arma::trans(rate_matrix->transition_rates(remap, migr_mat)), states, duration)
  {
    // loglikelihood of counts of coalescent events within window, with reduced emission and starting states
    // update states (probabilities that uncoalesced lineages are in XXX configuration)

    if(!(states.n_rows == num_trans + num_coal && states.n_cols == num_start)) Rcpp::stop("invalid states");
    if(!(n.n_elem == num_start)) Rcpp::stop("invalid n");
    if(!(y.n_rows == num_emiss && y.n_cols == num_start && arma::all(n - arma::trans(arma::sum(y,0)) >= 0))) Rcpp::stop("invalid y");

    states = transition.result;
    tran_probs = states.submat(arma::span(0, num_trans-1), arma::span::all);
    coal_probs = states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all);
    emiss_probs = arma::mat(num_emiss, states.n_cols, arma::fill::zeros);

    for (unsigned i=0; i<num_coal; ++i)
    {
      for(unsigned j=0; j<num_start; ++j)
      {
        emiss_probs.at(c_map_e(i,j),j) += coal_probs.at(i,j);
      }
    }
    arma::rowvec uncoal_probs = arma::ones<arma::rowvec>(num_start) - arma::sum(emiss_probs, 0);

    // if requested simulate new data
    if (simulate)
    {
      arma::imat y_sim (num_emiss+1, num_start);
      arma::mat probs = arma::join_vert(emiss_probs, uncoal_probs);
      for (unsigned i=0; i<num_start; ++i)
      {
        R::rmultinom(n.at(i), probs.colptr(i), num_emiss+1, y_sim.colptr(i));
      }
      y = arma::conv_to<arma::umat>::from(y_sim.head_rows(num_emiss));
    }

    // mask coalescent events that are not allowed for a given starting configuration
    for (unsigned i=0; i<num_start; ++i)
    {
      unsigned drop = 0;
      for (unsigned j=0; j<num_emiss; ++j)
      {
        if (s_mask_e.at(j,i))
        {
          drop += y.at(j,i);
          y.at(j,i) = 0;
        }
      }
      n.at(i) -= drop;
      // TODO warning about modifying input data?
    }
    y_in = y; // modify input y, as is done for input n

    // multinomial log-likelihood (missing data should be 0 in both n and y)
    z = n - arma::trans(arma::sum(y,0)); //member
    double lp = 0;
    double constant = 0;
    for (unsigned i=0; i<num_start; ++i)
    {
      unsigned nlp = 0;
      double psum = 0.;
      nlp += z.at(i);
      psum += uncoal_probs.at(i);
      lp += z.at(i)*log(uncoal_probs.at(i));
      if (normalizing_constant) constant += -1. * ::Rf_lgammafn(1+z.at(i));
      for (unsigned j=0; j<num_emiss; ++j)
      {
        if (s_mask_e.at(j,i) == 0)
        {
          nlp += y.at(j,i);
          psum += emiss_probs.at(j,i);
          lp += y.at(j,i)*log(emiss_probs.at(j,i));
          if (normalizing_constant) constant += -1. * ::Rf_lgammafn(1+y.at(j,i));
        }
      }
      //do I need to normalize to 1? eg will numerical errors cause total 
      if (std::fabs(1. - psum) >= sqrt(arma::datum::eps)) Rcpp::Rcout << "fixme:state probability sum " << psum << "\n";
      if (normalizing_constant) constant += 1. * ::Rf_lgammafn(1+nlp);
    }
    lp += constant;

    // condition on not coalescing
    //tran_probs.each_col([](arma::vec& x){ x /= arma::accu(x) ;});
    tran_probs.each_row([&](arma::rowvec& x){ x /= uncoal_probs ;});
    states.submat(arma::span(0, num_trans-1), arma::span::all) = tran_probs;
    states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all).zeros();

    // remaining lineages are uncoalesced lineages
    n = z;

    loglikelihood = lp; //member
  }

  arma::mat reverse_differentiate (arma::mat& d_states) 
  {
    // returns gradient of rate matrix with regard to multinomial loglikelihood
    // updates gradient of state vectors with regard to multinomial loglikelihood

    if(!(d_states.n_rows == num_trans + num_coal && d_states.n_cols == num_start)) Rcpp::stop("invalid d_states");

    arma::rowvec uncoal_probs = arma::ones<arma::rowvec>(num_start) - arma::sum(emiss_probs, 0);

    // conditional tran_probs --> unconditional tran_probs
    arma::mat d_tran_probs = d_states.submat(arma::span(0, num_trans-1), arma::span::all); 
    arma::rowvec d_uncoal_probs = -arma::sum(tran_probs % d_tran_probs, 0) / uncoal_probs;
    for (unsigned i=0; i<d_tran_probs.n_rows; ++i)
    {
      for (unsigned j=0; j<d_tran_probs.n_cols; ++j)
      {
        d_tran_probs.at(i,j) = d_tran_probs.at(i,j) / uncoal_probs.at(j);
      }
    }

    // loglikelihood --> emission probabilities
    arma::mat d_emiss_probs = arma::zeros(arma::size(emiss_probs));
    for (unsigned i=0; i<num_start; ++i)
    {
      for (unsigned j=0; j<num_emiss; ++j)
      {
        if (s_mask_e.at(j,i) == 0)
        {
          d_emiss_probs.at(j,i) += double(y.at(j,i))/emiss_probs.at(j,i);
        }
      }
      d_uncoal_probs.at(i) += double(z.at(i))/uncoal_probs.at(i);
    }
    d_emiss_probs.each_row([&](arma::rowvec& x) { x -= d_uncoal_probs; });

    // emission states --> coalescent states
    arma::mat d_coal_probs = arma::zeros(num_coal, num_start); 
    for (int i=int(num_coal)-1; i>=0; --i)
    {
      for (int j=int(num_start)-1; j>=0; --j)
      {
        d_coal_probs.at(i,j) += d_emiss_probs.at(c_map_e(i,j),j);
      }
    }

    // transition/coalescent states --> all states
    d_states.submat(arma::span(0, num_trans-1), arma::span::all) = d_tran_probs;
    d_states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all) = d_coal_probs;

    // current timepoint --> previous timepoint
    arma::mat d_states_updated, d_rates;
    transition.reverse_differentiate(d_rates, d_states_updated, d_states);
    d_states = d_states_updated;

    return rate_matrix->reverse_differentiate(arma::trans(d_rates), arma::trans(transition()), remap, migr_mat);
  }
};

// [[Rcpp::export()]]
Rcpp::List test_coalescent_epoch (arma::mat states, arma::uvec n, arma::umat y, arma::uvec remap, arma::mat migr_mat, arma::mat gradient)
{
  arma::mat new_states = states;
  transition_rate_matrix rate_matrix (migr_mat.n_rows, 3, true);
  double duration = 1.3;
  coalescent_epoch epoch (new_states, n, y, remap, migr_mat, duration, &rate_matrix);
  double ll = epoch.loglikelihood;
  arma::mat grad = epoch.reverse_differentiate (gradient);

  return Rcpp::List::create(
      Rcpp::_["states"] = new_states, 
      Rcpp::_["n"] = n, 
      Rcpp::_["grad"] = grad,
      Rcpp::_["gradient"] = gradient,
      Rcpp::_["rates"] = arma::trans(epoch.transition()), 
      Rcpp::_["ll"] = ll
      );
}

struct decoder 
{
  transition_rate_matrix trans_mat;

  decoder (const unsigned num_pop, const unsigned num_ind, const bool approx) : trans_mat(num_pop,num_ind,approx) {}

  arma::umat transitory_states (void)
  {
    return trans_mat.t_states;
  }

  arma::umat coalescent_states (void)
  {
    return trans_mat.c_states;
  }

  arma::umat emission_classes (void)
  {
    // rows are starting states, columns are populations, values are haplotype counts
    // these classes correspond to columns in data "y" / elements of "z"
    // used in loglikelihood()

    return trans_mat.s_states;
  }

  arma::umat emission_states (void)
  {
    // rows are emission types (types of pairwise coalescences), columns are first/second of pair, values are population indices
    // these states correspond to rows in data "y" used in loglikelihood()

    return trans_mat.e_states;
  }

  arma::umat map_coalescent_states_to_emission_states (void)
  {
    return trans_mat.c_map_e;
  }

  arma::uvec map_emission_classes_to_transitory_states (void)
  {
    return trans_mat.s_map_t;
  }

  std::string symbolic_transition_rates (const arma::uvec& remap)
  {
    return trans_mat.symbolic_transition_rates(remap);
  }

  arma::mat transition_probabilities (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    unsigned dim = trans_mat.num_trans + trans_mat.num_coal;
    arma::mat eye = arma::eye(dim, dim);
    matrix_exponential_multiply transition (arma::trans(trans_mat.transition_rates(remap, migr_mat)), eye, duration);
    return transition.result;
  }

  arma::mat transition_rates (const arma::uvec& remap, const arma::mat& migr_mat)
  {
    return trans_mat.transition_rates(remap, migr_mat);
  }

  arma::cube migration_function (const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(migr_mat.n_slices == T && migr_mat.n_rows == migr_mat.n_cols && arma::all(arma::vectorise(migr_mat) > 0.))) Rcpp::stop("invalid migr_mat");

    arma::cube migr_mat_copy = migr_mat;
    arma::cube migr_fun = arma::zeros<arma::cube>(migr_mat.n_rows, migr_mat.n_cols, T + 1);
    migr_fun.slice(0).eye();
    for (unsigned t=0; t<T; ++t)
    {
      for (unsigned i=0; i<migr_mat.n_rows; ++i)
      {
        migr_mat_copy.slice(t).at(i,i) = 0.;
        migr_mat_copy.slice(t).at(i,i) = -arma::accu(migr_mat_copy.slice(t).row(i));
      }
      matrix_exponential_multiply transition (arma::trans(migr_mat_copy.slice(t)), migr_fun.slice(t), duration.at(t));
      migr_fun.slice(t+1) = transition.result;
    }

    return migr_fun;
  }

  arma::cube migration_operator (const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(migr_mat.n_slices == T && migr_mat.n_rows == migr_mat.n_cols && arma::all(arma::vectorise(migr_mat) > 0.))) Rcpp::stop("invalid migr_mat");

    arma::cube migr_mat_copy = migr_mat;
    arma::cube migr_fun = arma::zeros<arma::cube>(migr_mat.n_rows, migr_mat.n_cols, T);
    for (unsigned t=0; t<T; ++t)
    {
      for (unsigned i=0; i<migr_mat.n_rows; ++i)
      {
        migr_mat_copy.slice(t).at(i,i) = 0.;
        migr_mat_copy.slice(t).at(i,i) = -arma::accu(migr_mat_copy.slice(t).row(i));
      }
      matrix_exponential_multiply transition (arma::trans(migr_mat_copy.slice(t)), arma::eye(migr_mat.n_rows, migr_mat.n_cols), duration.at(t));
      migr_fun.slice(t) = transition.result;
    }

    return migr_fun;
  }

  Rcpp::List migration_function_penalty (const arma::mat& states, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration, const arma::vec& penalty, const unsigned order)
  {
    unsigned T = duration.n_elem;

    if(!(migr_mat.n_slices == T && migr_mat.n_rows == migr_mat.n_cols && arma::all(arma::vectorise(migr_mat) >= 0.))) Rcpp::stop("invalid migr_mat");
    if(!(states.n_rows == migr_mat.n_cols && states.n_rows == states.n_cols && arma::all(arma::vectorise(states) >= 0.))) Rcpp::stop("invalid states");
    if(!(arma::all(penalty >= 0.) && penalty.n_elem == 5)) Rcpp::stop("invalid penalty");

    std::vector<penalized_migration_function> lin_hist;

    arma::cube migr_fun = arma::zeros<arma::cube>(migr_mat.n_rows, migr_mat.n_cols, T + 1);
    arma::mat states_copy = states; 
    arma::mat diff;

    double lp = 0.;
    migr_fun.slice(0) = states_copy;
    for (unsigned t=0; t<T; ++t)
    {
      lin_hist.emplace_back(states_copy, remap.col(t), migr_mat.slice(t), duration.at(t), penalty.at(3)); //last argument is another type of dirichlet sparsity
      migr_fun.slice(t+1) = lin_hist[t].migr_fun;
      lp += lin_hist[t].loglikelihood;
    }

    arma::cube d_migr_mat (arma::size(migr_mat), arma::fill::zeros);
    arma::cube d_migr_fun (arma::size(migr_fun), arma::fill::zeros);

    // gradient wrt to smoothness penalty on log migration rates
    if (penalty.at(0) > 0.)
    {
      diff = arma::eye(T, T);
      for (unsigned i=0; i<order; ++i)
      {
        arma::mat d (diff.n_rows - 1, diff.n_rows, arma::fill::zeros);
        d.diag().fill(-1.0);
        for(unsigned j=0; j<d.n_rows; ++j)
        {
          d.at(j,j+1) = 1.;
        }
        diff = d * diff;
      }
      for (unsigned i=0; i<migr_fun.n_rows; ++i)
      {
        for (unsigned j=0; j<migr_fun.n_cols; ++j)
        {
          arma::vec diff_ij = diff * arma::log10(arma::vectorise(migr_mat.tube(i,j)));
          lp += arma::accu(-0.5 * arma::pow(diff_ij, 2) * std::pow(penalty.at(0), 2));
          diff_ij = -1. * diff_ij * std::pow(penalty.at(0), 2);
          d_migr_mat.tube(i,j) = diff.t() * diff_ij; 
          d_migr_mat.tube(i,j) /= (migr_mat.tube(i,j) * log(10));
        }
      }
    }

    // gradient wrt to smoothness penalty on migration function
    if (penalty.at(1) > 0.)
    {
      diff = arma::eye(T+1, T+1);
      for (unsigned i=0; i<order; ++i)
      {
        arma::mat d (diff.n_rows - 1, diff.n_rows, arma::fill::zeros);
        d.diag().fill(-1.0);
        for(unsigned j=0; j<d.n_rows; ++j)
        {
          d.at(j,j+1) = 1.;
        }
        diff = d * diff;
      }
      for (unsigned i=0; i<migr_fun.n_rows; ++i)
      {
        for (unsigned j=0; j<migr_fun.n_cols; ++j)
        {
          arma::vec diff_ij = diff * arma::vectorise(migr_fun.tube(i,j));
          lp += arma::accu(-0.5 * arma::pow(diff_ij, 2) * std::pow(penalty.at(1), 2));
          diff_ij = -1. * diff_ij * std::pow(penalty.at(1), 2);
          d_migr_fun.tube(i,j) = diff.t() * diff_ij;
        }
      }
    }

    // penalty on sparsity of cumulative ancestry curves
    if (penalty.at(2) > 0.)
    {
      arma::mat migr_fun_total = arma::zeros(migr_mat.n_rows, migr_mat.n_cols),
                d_migr_fun_total = arma::zeros(migr_mat.n_rows, migr_mat.n_cols);
      for (unsigned t=0; t<T; ++t)
      {
        migr_fun_total += duration.at(t) * migr_fun.slice(t);
      }
      for (unsigned i=0; i<migr_fun_total.n_cols; ++i)
      {
        double colsum = arma::accu(migr_fun_total.col(i));
        arma::vec log_states = penalty.at(2) * arma::log(migr_fun_total.col(i)/colsum);
        double log_max = log_states.max();
        double log_sum = log_max + log(arma::accu(arma::exp(log_states - log_max)));
        lp += log_sum;
        arma::vec tmp = arma::exp(log(penalty.at(2)) + (penalty.at(2) - 1.) * arma::log(migr_fun_total.col(i)/colsum) - log_sum);
        double d_colsum = -arma::accu(tmp % migr_fun_total.col(i))/std::pow(colsum,2);
        d_migr_fun_total.col(i) = tmp / colsum + d_colsum;
      }
      for (int t=T-1; t>=0; --t)
      {
        d_migr_fun.slice(t) += duration.at(t) * d_migr_fun_total;
      }
    }

    // penalty on absolute size of migration rates
    if (penalty.at(4) > 0.)
    {
      // log(exp(-penalty * rate) * C) = -penalty * rate
      // d/log10 rate = -penalty * log(10) * 10^{log10 rate}
      for (unsigned i=0; i<migr_mat.n_rows; ++i)
      {
        for (unsigned j=0; j<migr_mat.n_rows; ++j)
        {
          if (i != j)
          {
            for (unsigned t=0; t<T; ++t)
            {
              lp += -penalty.at(4) * migr_mat.at(i,j,t);
              d_migr_mat.at(i,j,t) += -penalty.at(4);
            }
          }
        }
      }
    }

    // backpropagate gradient
    states_copy.zeros();
    for (int t=T-1; t>=0; --t)
    {
      d_migr_mat.slice(t) += lin_hist[t].reverse_differentiate(states_copy, d_migr_fun.slice(t+1));
    }

    return Rcpp::List::create(
        Rcpp::_["prediction"] = migr_fun, 
        Rcpp::_["deviance"] = -2 * lp,
        Rcpp::_["gradient"] = -2 * d_migr_mat,
        Rcpp::_["states"] = -2 * states_copy
        );
  }

  Rcpp::List debug (const arma::mat& states, const arma::uvec& n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(y.n_slices == T)) Rcpp::stop("invalid y");
    if(!(remap.n_cols == T)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == T)) Rcpp::stop("invalid migr_mat");
    if(!(duration.n_elem == T)) Rcpp::stop("invalid duration");

    std::vector<coalescent_epoch> epochs;
    epochs.reserve(T);

    arma::mat states_copy = states;
    arma::uvec n_copy = n;
    arma::cube emiss_probs (trans_mat.num_emiss + 1, trans_mat.num_start, T);
    arma::ucube y_corr (trans_mat.num_emiss + 1, trans_mat.num_start, T);
    arma::cube state_vectors (states.n_rows, states.n_cols, T+1);
    arma::vec loglik (T);

    state_vectors.slice(0) = states;

    double lp = 0.;
    for (unsigned t=0; t<T; ++t)
    {
      arma::umat y_copy = y.slice(t);
      epochs.emplace_back(states_copy, n_copy, y_copy, remap.col(t), migr_mat.slice(t), duration.at(t), &trans_mat);
      loglik.at(t) = epochs[t].loglikelihood;
      state_vectors.slice(t) = states_copy;
      lp += epochs[t].loglikelihood;
      y_corr.slice(t) = arma::join_vert(y_copy, arma::trans(epochs[t].z));
      emiss_probs.slice(t) = 
        arma::join_vert(epochs[t].emiss_probs, arma::ones<arma::rowvec>(trans_mat.num_start) - arma::sum(epochs[t].emiss_probs, 0));
    }

    arma::cube gradient (arma::size(migr_mat));
    states_copy.zeros();

    for (int t=T-1; t>=0; --t)
    {
      gradient.slice(t) = epochs[t].reverse_differentiate(states_copy);
    }

    return Rcpp::List::create(
        Rcpp::_["deviance"] = -2. * lp,
        Rcpp::_["gradient"] = -2. * gradient,
        Rcpp::_["states"] = -2. * states_copy,
        Rcpp::_["state_vectors"] = state_vectors,
        Rcpp::_["loglikelihood"] = -loglik,
        Rcpp::_["emission"] = y_corr,
        Rcpp::_["predicted"] = emiss_probs
        );
  }

  Rcpp::List deviance (const arma::mat& states, const arma::uvec& n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(y.n_slices == T)) Rcpp::stop("invalid y");
    if(!(remap.n_cols == T)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == T)) Rcpp::stop("invalid migr_mat");
    if(!(duration.n_elem == T)) Rcpp::stop("invalid duration");

    std::vector<coalescent_epoch> epochs;
    epochs.reserve(T);

    arma::mat states_copy = states;
    arma::uvec n_copy = n;
    arma::cube emiss_probs (trans_mat.num_emiss + 1, trans_mat.num_start, T);
    arma::ucube y_corr (trans_mat.num_emiss + 1, trans_mat.num_start, T);

    double lp = 0.;
    for (unsigned t=0; t<T; ++t)
    {
      arma::umat y_copy = y.slice(t);
      epochs.emplace_back(states_copy, n_copy, y_copy, remap.col(t), migr_mat.slice(t), duration.at(t), &trans_mat);
      lp += epochs[t].loglikelihood;
      y_corr.slice(t) = arma::join_vert(y_copy, arma::trans(epochs[t].z));
      emiss_probs.slice(t) = 
        arma::join_vert(epochs[t].emiss_probs, arma::ones<arma::rowvec>(trans_mat.num_start) - arma::sum(epochs[t].emiss_probs, 0));
    }

    arma::cube gradient (arma::size(migr_mat));
    states_copy.zeros();

    for (int t=T-1; t>=0; --t)
    {
      gradient.slice(t) = epochs[t].reverse_differentiate(states_copy);
    }

    return Rcpp::List::create(
        Rcpp::_["deviance"] = -2. * lp,
        Rcpp::_["gradient"] = -2. * gradient,
        Rcpp::_["states"] = -2. * states_copy,
        Rcpp::_["emission"] = y_corr,
        Rcpp::_["predicted"] = emiss_probs
        );
  }

  arma::ucube simulate (const arma::mat& states, const arma::uvec& n, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(remap.n_cols == T)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == T)) Rcpp::stop("invalid migr_mat");
    if(!(duration.n_elem == T)) Rcpp::stop("invalid duration");

    arma::mat states_copy = states;
    arma::uvec n_copy = n;
    arma::ucube y_sim (trans_mat.num_emiss, trans_mat.num_start, T);

    for (unsigned t=0; t<T; ++t)
    {
      arma::umat y = arma::zeros<arma::umat>(trans_mat.num_emiss, trans_mat.num_start);
      coalescent_epoch epoch (states_copy, n_copy, y, remap.col(t), migr_mat.slice(t), duration.at(t), &trans_mat, true);
      y_sim.slice(t) = y;
    }

    return y_sim;
  }

  arma::mat initial_states (void)
  {
    arma::mat states = arma::eye<arma::mat>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans);
    return states.cols(trans_mat.s_map_t);
  }
};

struct transition_matrix
{
  unsigned P, I;
  arma::umat t_states, c_states, e_states, s_states, c_map_e;
  arma::uvec s_map_t;
  unsigned num_trans, num_coal, num_start, num_emiss;

  transition_matrix (const unsigned P, const unsigned I) : P(P), I(I)
  {
    // enumerate transitory states (e.g. haplotype configurations across populations)
    arma::uvec pop_index = arma::regspace<arma::uvec>(0, P-1);
    t_states = pop_index;
    for (unsigned i=1; i<I; ++i)
      t_states = arma::join_horiz(arma::repmat(t_states, pop_index.n_elem, 1), arma::repelem(pop_index, t_states.n_rows, 1));

    // find nonredundant starting configurations
    s_map_t = arma::regspace<arma::uvec>(0, t_states.n_rows-1);
    std::sort(s_map_t.begin(), s_map_t.end(), [&](unsigned i, unsigned j){
              arma::urowvec a = arma::sort(t_states.row(i)),
                            b = arma::sort(t_states.row(j));
              for(unsigned k=0; k<I-1; ++k) 
                if(a.at(k) != b.at(k))
                   return a.at(k) < b.at(k);
              return a.at(I-1) < b.at(I-1);
    });
    auto end = std::unique(s_map_t.begin(), s_map_t.end(), [&](unsigned i, unsigned j){
                arma::irowvec a = arma::conv_to<arma::irowvec>::from(arma::sort(t_states.row(i))),
                              b = arma::conv_to<arma::irowvec>::from(arma::sort(t_states.row(j)));
                if (arma::any(a - b))
                  return false;
                return true;
    });
    s_map_t = s_map_t.head(end - s_map_t.begin());
    s_states = arma::umat(s_map_t.n_elem, P, arma::fill::zeros);
    for (unsigned i=0; i<s_map_t.n_elem; ++i)
      for (unsigned j=0; j<I; ++j)
        s_states.at(i,t_states.at(s_map_t.at(i),j)) += 1;

    // enumerate possible pairwise coalescence events
    c_states = arma::zeros<arma::umat>(0, 2);
    for (unsigned i=0; i<I; ++i)
      for (unsigned j=i; j<I; ++j)
        if (i != j)
          c_states = arma::join_vert(c_states, arma::urowvec({i,j}));

    // find nonredundant coalescent events for a given starting state
    // e.g. conditional on starting state, divide coal events between indexed lineages into within/cross population
    e_states = arma::zeros<arma::umat>(0, 2);
    for (unsigned i=0; i<P; ++i)
      for (unsigned j=i; j<P; ++j)
        e_states = arma::join_vert(e_states, arma::urowvec({i,j}));
    c_map_e = arma::zeros<arma::umat>(c_states.n_rows, s_map_t.n_elem);
    for(unsigned i=0; i<s_map_t.n_elem; ++i)
      for(unsigned j=0; j<c_states.n_rows; ++j)
      {
        unsigned p = t_states.at(s_map_t.at(i),c_states.at(j,0)), 
                 q = t_states.at(s_map_t.at(i),c_states.at(j,1));
        arma::uvec idx = arma::find(e_states.col(0) == std::min(p,q) && e_states.col(1) == std::max(p,q), 1);
        if (idx.n_elem)
          c_map_e.at(j,i) = idx.at(0);
        else
          Rcpp::warning("coalescent event does not map to emission");
      }

    num_trans = t_states.n_rows;
    num_coal  = c_states.n_rows;
    num_start = s_states.n_rows;
    num_emiss = e_states.n_rows;
  }

  std::string transition_rates_symbolic (const arma::uvec& remap)
  {

    if (!(remap.n_elem == P && arma::max(remap) < P)) Rcpp::stop("invalid remap");

    std::string out = "";
    std::string outstring;

    // transition rates to transient states
    for (unsigned i=0; i<num_trans; ++i)
      for (unsigned j=0; j<num_trans; ++j)
      {
        if (i != j)
        {
          unsigned moves = 0;
          outstring  = "M[" + std::to_string(i+1) + "," + std::to_string(j+1) + "] = 1";
          for (unsigned k=0; k<I; ++k)
          {
            unsigned p1 = remap.at(t_states.at(i,k)), // origin, backwards in time; destination, forwards in time
                     p2 = remap.at(t_states.at(j,k)); // destination, backwards in time; origin, forwards in time
            if (p1 != p2)
            {
              moves += 1;
              outstring += "*m[" + std::to_string(p1+1) + "," + std::to_string(p2+1) + "]";
            }
          }
          if (moves == 0) 
            outstring = "M[" + std::to_string(i+1) + "," + std::to_string(j+1) + "] = 0";
          outstring += ";\n";
          out += outstring;
        }
      }

    // transition rates to coalesced states
    for (unsigned i=0; i<num_trans; ++i)
    {
      for (unsigned j=0; j<num_coal; ++j)
      {
        const unsigned p1 = remap.at(t_states.at(i,c_states.at(j,0))),
                       p2 = remap.at(t_states.at(i,c_states.at(j,1)));
        outstring = "M[" + std::to_string(i+1) + "," + std::to_string(num_trans+j+1) + "]";
        outstring += p1==p2 ? 
          " = 1/m[" + std::to_string(p1+1) + "," + std::to_string(p1+1) + "];\n" :
          " = 0;\n";
        out += outstring;
      }
    }

    // total exit rates
    for (unsigned i=0; i<num_trans+num_coal; ++i)
    {
      outstring = "M[" + std::to_string(i+1) + "," + std::to_string(i+1) + "] = ";
      for (unsigned j=0; j<num_trans+num_coal; ++j)
        if (i != j)
          outstring += "-M[" + std::to_string(i+1) + "," + std::to_string(j+1) + "]";
      outstring += ";\n";
      out += outstring;
    }

   return out;
  }

  arma::mat transition_rates (const arma::uvec& remap, arma::mat migr_mat, const double duration)
  {
    // Transition rate matrix for migration process up until first coalescent event
    // with populations remapped to a population vector of lower order. This results in a transition
    // matrix with redundant states. 

    arma::vec ne = migr_mat.diag();
    migr_mat.diag().ones();

    if (!(remap.n_elem == P && arma::max(remap) < P)) Rcpp::stop("invalid remap");
    if (!(migr_mat.n_rows == P && migr_mat.n_cols == P && arma::all(arma::vectorise(migr_mat) >= 0.))) Rcpp::stop("invalid migr_mat");
    if (!(arma::all(ne >= 0.))) Rcpp::stop("invalid ne");
    if (!(duration >= 0.)) Rcpp::stop("invalid duration");

    arma::mat rates (num_trans+num_coal, num_trans+num_coal, arma::fill::zeros);

    // transition rates to transient states
    for (unsigned i=0; i<num_trans; ++i)
      for (unsigned j=0; j<num_trans; ++j)
      {
        unsigned moves = 0;
        rates.at(i,j) = 1.;
        for (unsigned k=0; k<I; ++k)
        {
          unsigned p1 = remap.at(t_states.at(i,k)), // origin, backwards in time; destination, forwards in time
                   p2 = remap.at(t_states.at(j,k)); // destination, backwards in time; origin, forwards in time
          if (p1 != p2)
          {
            moves += 1;
            rates.at(i,j) *= migr_mat.at(p1,p2);
            // dividing by ne here means that migration is in general asymmetric, but puts all parameters on the same scale (number haploids)
            // might want to reconsider this parameterization if models with symmetric migration are desired
          }
        }
        if (moves == 0) rates.at(i,j) = 0.;
      }

    // transition rates to coalesced states
    for (unsigned i=0; i<num_trans; ++i)
      for (unsigned j=0; j<num_coal; ++j)
      {
        const unsigned p1 = remap.at(t_states.at(i,c_states.at(j,0))),
                       p2 = remap.at(t_states.at(i,c_states.at(j,1)));
        rates.at(i,num_trans+j) = p1==p2 ? 1./ne[p1] : 0.;
      }

    // total exit rates
    for (unsigned i=0; i<rates.n_rows; ++i)
    {
      rates.at(i,i) = 0.;
      for (unsigned j=0; j<rates.n_cols; ++j)
        if (i != j)
          rates.at(i,i) -= rates.at(i,j);
    }

   return rates;
  }

  arma::mat transition_probabilities (const arma::mat& migr_mat, const double duration)
  {
    arma::uvec remap = arma::regspace<arma::uvec>(0, P-1);
    return transition_probabilities(remap, migr_mat, duration);
  }

  arma::mat transition_probabilities (const arma::uvec& remap, arma::mat migr_mat, const double duration)
  {
    // Transition probability matrix for migration process up until first coalescent event
    // with populations remapped to a population vector of lower order. This results in a transition
    // matrix with redundant states. 

    arma::vec ne = migr_mat.diag();
    migr_mat.diag().ones();

    if (!(remap.n_elem == P && arma::max(remap) < P)) Rcpp::stop("invalid remap");
    if (!(migr_mat.n_rows == P && migr_mat.n_cols == P && arma::all(arma::vectorise(migr_mat) >= 0.))) Rcpp::stop("invalid migr_mat");
    if (!(arma::all(ne >= 0.))) Rcpp::stop("invalid ne");
    if (!(duration >= 0.)) Rcpp::stop("invalid duration");

    arma::mat rates (num_trans+num_coal, num_trans+num_coal, arma::fill::zeros);

    // transition rates to transient states
    for (unsigned i=0; i<num_trans; ++i)
      for (unsigned j=0; j<num_trans; ++j)
      {
        unsigned moves = 0;
        rates.at(i,j) = 1.;
        for (unsigned k=0; k<I; ++k)
        {
          unsigned p1 = remap.at(t_states.at(i,k)),
                   p2 = remap.at(t_states.at(j,k));
          if (p1 != p2)
          {
            moves += 1;
            rates.at(i,j) *= migr_mat.at(p1,p2);//TODO ... should it be p2 first if parameterizing in terms of forward migration rates
          }
        }
        if (moves == 0) rates.at(i,j) = 0.;
      }

    // transition rates to coalesced states
    for (unsigned i=0; i<num_trans; ++i)
      for (unsigned j=0; j<num_coal; ++j)
      {
        const unsigned p1 = remap.at(t_states.at(i,c_states.at(j,0))),
                       p2 = remap.at(t_states.at(i,c_states.at(j,1)));
        rates.at(i,num_trans+j) = p1==p2 ? 1./ne[p1] : 0.;
      }

    // total exit rates
    for (unsigned i=0; i<rates.n_rows; ++i)
    {
      rates.at(i,i) = 0.;
      for (unsigned j=0; j<rates.n_cols; ++j)
        if (i != j)
          rates.at(i,i) -= rates.at(i,j);
    }

    arma::mat probs (arma::size(rates), arma::fill::zeros);
    if (!arma::expmat(probs, duration * rates) ||
        probs.has_inf() || 
        probs.has_nan() || 
        arma::any(arma::vectorise(probs) < 0.))
    {
      //Rcpp::Rcout << "arma::expmat failed, using leaky expm::expm" << std::endl;
      probs = expm(duration * rates);
    }
    
    //arma::mat probs = expm(duration * rates); //this has memory leak somehow
    //if (probs.has_inf() || probs.has_nan() || arma::any(arma::vectorise(probs) < 0.))
    //{
    //  Rcpp::Rcout << "reverting to an even more dubious way to calculate matrix exponential" << std::endl;
    //  //migr_mat.print("migr_mat");
    //  //ne.print("ne");
    //  probs = arma::expmat(duration * rates);
    //  //probs.print("em");
    //}

    return probs;
  }

  double loglikelihood (arma::mat& states, arma::uvec& n, const arma::umat& y, const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return loglikelihood (states, n, y, transition_probabilities(remap, migr_mat, duration));
  }

  double loglikelihood (arma::mat& states, arma::uvec& n, const arma::umat& y, const arma::mat& transition)
  {
    // loglikelihood of counts of coalescent events within window, with reduced emission and starting states
    // update states (probabilities that uncoalesced lineages are in XXX configuration)

    if(!(states.n_rows == num_trans + num_coal && states.n_cols == num_start)) Rcpp::stop("invalid states");
    if(!(n.n_elem == num_start)) Rcpp::stop("invalid n");
    if(!(y.n_rows == num_emiss && y.n_cols == num_start && arma::all(n - arma::trans(arma::sum(y,0)) >= 0))) Rcpp::stop("invalid y");
    if(!(transition.n_rows == num_trans+num_coal && transition.n_cols == transition.n_rows && arma::all(arma::vectorise(transition) >= 0.))) Rcpp::stop("invalid transition");

    states = arma::trans(transition) * states;
    arma::mat tran_probs = states.submat(arma::span(0, num_trans-1), arma::span::all),
              coal_probs = states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all);

    arma::mat emiss_probs = arma::mat(num_emiss, states.n_cols, arma::fill::zeros);
    for (unsigned i=0; i<num_coal; ++i)
      for(unsigned j=0; j<num_start; ++j)
        emiss_probs.at(c_map_e(i,j),j) += coal_probs.at(i,j);

    arma::rowvec uncoal_probs = arma::ones<arma::rowvec>(num_start) - arma::sum(emiss_probs, 0);

    // multinomial log-likelihood (missing data should be 0 in both n and y)
    arma::uvec z = n - arma::trans(arma::sum(y,0));
    double lp = 0;
    for (unsigned i=0; i<num_start; ++i)
    {
      unsigned nlp = 0;
      if(uncoal_probs.at(i) > 0.0 || z.at(i) > 0)
      {
        nlp += z.at(i);
        lp  += z.at(i)*log(uncoal_probs.at(i)) - ::Rf_lgammafn(1+z.at(i));
      }
      for (unsigned j=0; j<num_emiss; ++j)
        if(emiss_probs.at(j,i) > 0.0 || y.at(j,i) > 0)
        {
          nlp += y.at(j,i);
          lp  += y.at(j,i)*log(emiss_probs.at(j,i)) - ::Rf_lgammafn(1+y.at(j,i));
        }
      lp += ::Rf_lgammafn(1+nlp);
    }

    // condition on not coalescing
    tran_probs.each_col([](arma::vec& x){ x /= arma::accu(x) ;});
    states.submat(arma::span(0, num_trans-1), arma::span::all) = tran_probs;
    states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all).zeros();

    // remaining lineages are uncoalesced lineages
    n = z;

    return lp;
  }

  arma::mat reverse_differentiate (arma::mat& d_states, arma::uvec& n, const arma::umat& y, matrix_exponential_multiply& transition,
                                   const arma::mat& tran_probs, const arma::rowvec& uncoal_probs, const arma::mat& emiss_probs)
  {
    // returns gradient of rate matrix with regard to multinomial loglikelihood
    // updates gradient of state vectors with regard to multinomial loglikelihood
    // updates tree count n BACKWARDS in time (e.g. it increases)

    if(!(d_states.n_rows == num_trans + num_coal && d_states.n_cols == num_start)) Rcpp::stop("invalid d_states");
    if(!(n.n_elem == num_start)) Rcpp::stop("invalid n");
    if(!(y.n_rows == num_emiss && y.n_cols == num_start && arma::all(n - arma::trans(arma::sum(y,0)) >= 0))) Rcpp::stop("invalid y");
    if(!(transition.dim == num_trans+num_coal)) Rcpp::stop("invalid transition");
    //checks for tran_probs etc if not saving in struct

    //save from forward:
    // tran_probs
    // uncoal_probs
    // emiss_probs

    arma::uvec z = n;
    n += arma::trans(arma::sum(y,0));

    // conditional tran_probs --> unconditional tran_probs
    arma::mat d_tran_probs = 
      d_states.submat(arma::span(0, num_trans-1), arma::span::all) % 
      (arma::ones(tran_probs.n_rows) * uncoal_probs) % (1 - tran_probs);

    // loglikelihood --> emission probabilities
    arma::mat d_emiss_probs = arma::zeros(arma::size(emiss_probs));
    arma::rowvec d_uncoal_probs = arma::zeros(arma::size(uncoal_probs));
    for (unsigned i=num_start-1; i>=0; --i)
    {
      for (unsigned j=num_emiss-1; j>=0; --j)
      {
        if(emiss_probs.at(j,i) > 0.0 || y.at(j,i) > 0)
        {
          d_emiss_probs.at(j,i) = y.at(j,i)/emiss_probs.at(j,i);
        }
      }
      if(uncoal_probs.at(i) > 0.0 || z.at(i) > 0)
      {
        d_uncoal_probs.at(i) = z.at(i)/uncoal_probs.at(i);
      }
    }
    d_emiss_probs += -arma::ones(d_emiss_probs.n_rows) * d_uncoal_probs.t(); 

    // emission states --> coalescent states
    arma::mat d_coal_probs = arma::zeros(num_coal, num_start); 
    for (unsigned i=num_coal-1; i>=0; --i)
    {
      for (unsigned j=num_start-1; j>=0; --j)
      {
        d_coal_probs.at(i,j) = d_emiss_probs.at(c_map_e(i,j),j);
      }
    }

    // transition/coalescent states --> all states
    d_states.submat(arma::span(0, num_trans-1), arma::span::all) = d_tran_probs;
    d_states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all) = d_coal_probs;

    // current timepoint --> previous timepoint
    arma::mat d_states_copy = d_states;
    arma::mat d_rates = arma::zeros(transition.dim, transition.dim);
    
    transition.reverse_differentiate(d_states, d_rates, d_states_copy); //d_states is persistant reference

    return d_rates;
  }

  arma::umat simulate (arma::mat& states, arma::uvec& n, const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return simulate (states, n, transition_probabilities(remap, migr_mat, duration));
  }
  
  arma::umat simulate (arma::mat& states, arma::uvec& n, const arma::mat& transition)
  {
    // simulate counts of types of coalescent events (which population-pair coalesces first) 
    // within window for each class (starting haplotype configuration)
    // the last row of output are the remaining (uncoalesced) lineages at end of window ("z" in loglikelihood)

    if(!(states.n_rows == num_trans + num_coal && states.n_cols == num_start)) Rcpp::stop("invalid states");
    if(!(n.n_elem == num_start)) Rcpp::stop("invalid n");
    if(!(transition.n_rows == num_trans+num_coal && transition.n_cols == transition.n_rows && arma::all(arma::vectorise(transition) >= 0.))) Rcpp::stop("invalid transition");

    states = arma::trans(transition) * states;
    arma::mat tran_probs = states.submat(arma::span(0, num_trans-1), arma::span::all),
              coal_probs = states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all);

    arma::mat emiss_probs = arma::mat(num_emiss, states.n_cols, arma::fill::zeros);
    for (unsigned i=0; i<num_coal; ++i)
      for(unsigned j=0; j<num_start; ++j)
        emiss_probs.at(c_map_e(i,j),j) += coal_probs.at(i,j);

    emiss_probs = arma::join_vert(emiss_probs, arma::ones<arma::rowvec>(num_start) - arma::sum(emiss_probs, 0));

    // multinomial simulation
    arma::imat y (num_emiss+1, num_start);
    for (unsigned i=0; i<num_start; ++i)
      R::rmultinom(n.at(i), emiss_probs.colptr(i), num_emiss+1, y.colptr(i));

    // condition on not coalescing
    tran_probs.each_col([](arma::vec& x){ x /= arma::accu(x) ;});
    states.submat(arma::span(0, num_trans-1), arma::span::all) = tran_probs;
    states.submat(arma::span(num_trans, num_trans+num_coal-1), arma::span::all).zeros();

    // remaining lineages are uncoalesced lineages
    n = arma::conv_to<arma::uvec>::from(arma::trans(y.row(y.n_rows-1)));

    return arma::conv_to<arma::umat>::from(y.head_rows(num_emiss));
  }
};

struct particle
{
  transition_matrix trans_mat;

  particle (const unsigned num_pop, const unsigned num_ind) : trans_mat(num_pop,num_ind) {}

  arma::umat transitory_states (void)
  {
    return trans_mat.t_states;
  }

  arma::umat coalescent_states (void)
  {
    return trans_mat.c_states;
  }

  arma::umat emission_classes (void)
  {
    // rows are starting states, columns are populations, values are haplotype counts
    // these classes correspond to columns in data "y" / elements of "z"
    // used in loglikelihood()

    return trans_mat.s_states;
  }

  arma::umat emission_states (void)
  {
    // rows are emission types (types of pairwise coalescences), columns are first/second of pair, values are population indices
    // these states correspond to rows in data "y" used in loglikelihood()

    return trans_mat.e_states;
  }

  arma::umat map_coalescent_states_to_emission_states (void)
  {
    return trans_mat.c_map_e;
  }

  arma::uvec map_emission_classes_to_transitory_states (void)
  {
    return trans_mat.s_map_t;
  }

  std::string symbolic_transition_rates (const arma::uvec& remap)
  {
    return trans_mat.transition_rates_symbolic(remap);
  }

  arma::mat transition_probabilities (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_probabilities(remap, migr_mat, duration);
  }

  arma::mat transition_rates (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_rates(remap, migr_mat, duration);
  }

  double partial (arma::mat& states, const arma::uvec& n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(y.n_slices == T)) Rcpp::stop("invalid y");
    if(!(remap.n_cols == T)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == T)) Rcpp::stop("invalid migr_mat");
    if(!(duration.n_elem == T)) Rcpp::stop("invalid duration");

    arma::uvec z = n;
    double lp = 0.;
    for (unsigned t=0; t<T; ++t)
      lp += trans_mat.loglikelihood(states, z, y.slice(t), remap.col(t), migr_mat.slice(t), duration.at(t));

    return lp;
  }

  double loglikelihood (const arma::uvec& n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    arma::mat states = arma::eye<arma::mat>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans);
    states = states.cols(trans_mat.s_map_t);
    return partial(states, n, y, remap, migr_mat, duration);
  }

  arma::ucube simulate (arma::uvec n, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    unsigned T = duration.n_elem;

    if(!(remap.n_cols == T)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == T)) Rcpp::stop("invalid migr_mat");
    if(!(duration.n_elem == T)) Rcpp::stop("invalid duration");

    double lp = 0.;
    arma::ucube y (trans_mat.num_emiss, trans_mat.num_start, T);
    arma::mat states = arma::eye<arma::mat>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans);
    states = states.cols(trans_mat.s_map_t);
    for (unsigned t=0; t<T; ++t)
      y.slice(t) = trans_mat.simulate(states, n, remap.col(t), migr_mat.slice(t), duration.at(t));

    return y;
  }
};

struct particulate
{
  // struct particle above is really inefficient for updating parameters at a single window
  // and calculating likelihood. this does that much more efficiently by caching the transition matrix
  // and state variables.

  transition_matrix trans_mat;

  arma::ucube _y;
  arma::umat _remap, _n;
  arma::vec  _duration, _loglik;
  arma::cube _parameters, _state, _transition;

  particulate (const unsigned num_pop, const unsigned num_ind)
  : trans_mat (num_pop, num_ind)
  {}

  arma::umat transitory_states (void)
  {
    return trans_mat.t_states;
  }

  arma::umat coalescent_states (void)
  {
    return trans_mat.c_states;
  }

  arma::umat emission_classes (void)
  {
    // rows are starting states, columns are populations, values are haplotype counts
    // these classes correspond to columns in data "y" / elements of "z"
    // used in loglikelihood()

    return trans_mat.s_states;
  }

  arma::umat emission_states (void)
  {
    // rows are emission types (types of pairwise coalescences), columns are first/second of pair, values are population indices
    // these states correspond to rows in data "y" used in loglikelihood()

    return trans_mat.e_states;
  }

  arma::mat transition_probabilities (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_probabilities(remap, migr_mat, duration);
  }

  arma::mat transition_rates (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_rates(remap, migr_mat, duration);
  }

  double init (arma::uvec n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration)
  {
    // calculate all transition matrices then loglikelihood
    // store n (starting lineages) for each window
    // also store _starting_ state for each window

    if(!(y.n_slices == duration.n_elem)) Rcpp::stop("invalid y");
    if(!(remap.n_cols == duration.n_elem)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == duration.n_elem)) Rcpp::stop("invalid migr_mat");

    _y = y;
    _remap = remap;
    _parameters = migr_mat;
    _duration = duration;
    _n = arma::zeros<arma::umat>(n.n_elem, duration.n_elem);
    _loglik = arma::zeros<arma::vec>(duration.n_elem);
    _transition = arma::zeros<arma::cube>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans+trans_mat.num_coal, duration.n_elem);
    _state = arma::zeros<arma::cube>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_start, duration.n_elem);

    arma::mat states = arma::eye<arma::mat>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans);
    states = states.cols(trans_mat.s_map_t);
    for (unsigned i=0; i<duration.n_elem; ++i)
    {
      _transition.slice(i) = trans_mat.transition_probabilities(_remap.col(i), _parameters.slice(i), _duration.at(i));
      _state.slice(i) = states;
      _n.col(i) = n;
      _loglik.at(i) = -trans_mat.loglikelihood(states, n, _y.slice(i), _transition.slice(i));
    }

    return arma::accu(_loglik);
  }

  double loglikelihood (void) const
  {
    return arma::accu(_loglik);
  }

  double update (const unsigned i, const arma::mat& parameters)
  {
    // calculate first transition matrix, use cached transition matrices
    // to recalculate _starting_ state, loglikelihood in downstream windows,

    if (i >= _duration.n_elem) Rcpp::stop("invalid i");

    _transition.slice(i) = trans_mat.transition_probabilities(_remap.col(i), parameters, _duration.at(i));
    _parameters.slice(i) = parameters;

    arma::mat states = _state.slice(i);
    arma::uvec n = _n.col(i);

    for (unsigned j=i; j<_duration.n_elem; ++j)
    {
      _state.slice(j) = states;
      _loglik.at(j) = -trans_mat.loglikelihood(states, n, _y.slice(j), _transition.slice(j));
    }
     
    // return loglikelihood from _entire_ sequence (including upstream)
    return arma::accu(_loglik);
  }

  arma::mat parameters (const unsigned i) const
  {
    return _parameters.slice(i);
  }

  arma::mat state (const unsigned i) const
  {
    return _state.slice(i);
  }

  arma::mat transition (const unsigned i) const
  {
    return _transition.slice(i);
  }
};

struct particoal
{
  // struct particle above is really inefficient for updating parameters at a single window
  // and calculating likelihood. this does that much more efficiently by caching the transition matrix
  // and state variables.

  transition_matrix trans_mat;

  arma::ucube _y, _mask;
  arma::umat _remap, _n;
  arma::vec  _duration, _loglik;
  arma::cube _parameters, _state, _transition;

  particoal (const unsigned num_pop, const unsigned num_ind)
  : trans_mat (num_pop, num_ind)
  {}

  arma::umat transitory_states (void)
  {
    // internal transitory states

    return trans_mat.t_states;
  }

  arma::umat coalescent_states (void)
  {
    // internal coalescent states

    return trans_mat.c_states;
  }

  arma::umat emission_classes (void)
  {
    // rows are starting states, columns are populations, values are haplotype counts
    // these classes correspond to columns in data "y" / elements of "z"
    // used in loglikelihood()

    return trans_mat.s_states;
  }

  arma::umat emission_states (void)
  {
    // rows are emission types (types of pairwise coalescences), columns are first/second of pair, values are population indices
    // these states correspond to rows in data "y" used in loglikelihood()

    return trans_mat.e_states;
  }

  arma::mat transition_probabilities (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_probabilities(remap, migr_mat, duration);
  }

  arma::mat transition_rates (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_rates(remap, migr_mat, duration);
  }

  struct internal_nlopt_data
  {
    // these should really be const
    unsigned i;
    particoal& p;

    // these are mutable
    int        err;
    unsigned   num_par;
    arma::vec  start, lb;
    arma::umat mapping;
    arma::mat  parameters, states;
    arma::uvec n;
    double     loglik;

    // settings
    bool verbose = true;
    bool smooth = false;

    static double filter (unsigned dim, const double* x, double* grad, void* f_data)
    {
      internal_nlopt_data* data = static_cast<internal_nlopt_data*>(f_data);

      unsigned i = data->i;

      data->n      = data->p._n.col(i);
      data->states = data->p._state.slice(i);

      if (arma::max(arma::vectorise(data->mapping)) > dim) Rcpp::stop("bad mask");

      data->parameters.zeros();
      for (unsigned j=0; j<data->mapping.n_elem; ++j)
        if (data->mapping.at(j) > 0)
          data->parameters.at(j) = x[data->mapping.at(j)-1];

      double lp = data->p.trans_mat.loglikelihood(data->states, data->n, data->p._y.slice(i), data->p._remap.col(i), data->parameters, data->p._duration.at(i));

      // if we're doing smoothing, we propagate the state vector forward and incorporate loglikelihood of future observations
      if (data->smooth)
        for (unsigned j=i+1; j<data->p._duration.n_elem; ++j)
          lp += data->p.trans_mat.loglikelihood(data->states, data->n, data->p._y.slice(j), data->p._transition.slice(j));

      return -lp;
    }

    internal_nlopt_data (particoal& p, const unsigned i, bool do_smooth, bool symmetric_migration) : i(i), p(p)
    {
      if (i >= p._duration.n_elem) Rcpp::stop("bad i");

      smooth = do_smooth;

      parameters = p._parameters.slice(i);

      // find set of active parameters and linear indices thereof
      arma::uvec active = arma::unique(p._remap.col(i));
      mapping = arma::zeros<arma::umat>(p._remap.n_rows, p._remap.n_rows);
      num_par = 0;
      if (symmetric_migration)
      {
        for (auto j : active)
          for (auto k : active)
          {
            if (k >= j)
            {
              num_par += 1;
              mapping.at(k,j) = mapping.at(j,k) = num_par;
            }
          }
      } else {
        for (auto j : active)
          for (auto k : active)
          {
            num_par += 1;
            mapping.at(k,j) = num_par;
          }
      }

      // lower bounds
      double min_ne = 1.;
      lb = arma::zeros<arma::vec>(num_par);
      unsigned l = 0;
      if (symmetric_migration)
      {
        for (auto j : active)
        {
          for (auto k : active)
          {
            if (k >= j)
              lb.at(l++) = k == j ? min_ne : 0.;
          }
        }
      } 
      else 
      {
        for (auto j : active)
        {
          for (auto k : active)
          {
            lb.at(l++) = k == j ? min_ne : 0.;
          }
        }
      }

      // get starting values
      start = std::vector<double>(num_par);
      for (unsigned j=0; j<parameters.n_elem; ++j)
      {
        if (mapping.at(j) > 0)
        {
          start[mapping.at(j)-1] = std::max(parameters.at(j), lb.at(mapping.at(j)-1));
        }
      }

      nlopt_opt optimizer = nlopt_create(NLOPT_LN_BOBYQA, num_par);
      nlopt_set_min_objective(optimizer, filter, this);
      nlopt_set_lower_bounds(optimizer, &(lb[0]));
      nlopt_set_xtol_rel(optimizer, 1e-4);
      err = nlopt_optimize(optimizer, &(start[0]), &loglik);

      if (err < 0)
      { // error; revert to initial state
        if (verbose) std::cout << "[nlopt] Failed with code: " << err << std::endl;
        parameters.print("Dumping parameters then reverting:");
        parameters = p._parameters.slice(i);
        for (unsigned j=0; j<parameters.n_elem; ++j)
          if (mapping.at(j) > 0)
            start[mapping.at(j)-1] = parameters.at(j);
      } else {
        for (unsigned j=0; j<parameters.n_elem; ++j)
          if (mapping.at(j) > 0)
            parameters.at(j) = start[mapping.at(j)-1];
      }

      // ensure that state/n are consistent with optimized parameters; and that loglik is just for the current observation
      smooth = false;
      loglik = filter(num_par, &(start[0]), nullptr, this);

      nlopt_destroy(optimizer);
    }
  };

  double estimate (arma::uvec n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration, const unsigned num_smooth)
  {
    // calculate all transition matrices then loglikelihood
    // store n (starting lineages) for each window
    // also store _starting_ state for each window

    bool symmetric_migration = true;

    if(!(y.n_slices == duration.n_elem)) Rcpp::stop("invalid y");
    if(!(remap.n_cols == duration.n_elem)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == duration.n_elem)) Rcpp::stop("invalid migr_mat");

    _y = y;
    _remap = remap;
    _parameters = migr_mat;
    _duration = duration;
    _n = arma::zeros<arma::umat>(n.n_elem, duration.n_elem);
    _loglik = arma::zeros<arma::vec>(duration.n_elem);
    _transition = arma::zeros<arma::cube>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans+trans_mat.num_coal, duration.n_elem);
    _state = arma::zeros<arma::cube>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_start, duration.n_elem);

    // filtering step
    arma::mat states = arma::eye<arma::mat>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans);
    states = states.cols(trans_mat.s_map_t);
    for (unsigned i=0; i<_duration.n_elem; ++i)
    {
      // inputs to optimization
      _state.slice(i) = states;
      _n.col(i) = n;

      internal_nlopt_data filter (*this, i, false, symmetric_migration);

      // update running variables
      n = filter.n;
      states = filter.states;

      // outputs of optimization
      _loglik.at(i) = filter.loglik;
      _parameters.slice(i) = filter.parameters;
      _transition.slice(i) = trans_mat.transition_probabilities(_remap.col(i), _parameters.slice(i), _duration.at(i));
    }
    std::cout << "Deviance after filtering: " << 2.*arma::accu(_loglik) << std::endl;

    // smoothing steps
    double old_loglik;
    for (unsigned s=0; s<num_smooth; ++s)
    {
      for (unsigned i=_duration.n_elem; i>0; --i)
      {
        internal_nlopt_data smooth (*this, i-1, true, symmetric_migration);

        // outputs of optimization
        _loglik.at(i-1)        = smooth.loglik;
        _parameters.slice(i-1) = smooth.parameters;
        _transition.slice(i-1) = trans_mat.transition_probabilities(_remap.col(i-1), _parameters.slice(i-1), _duration.at(i-1));

        // propagate the new states forward and update loglikelihood of future timepoints
        states = smooth.states;
        n      = smooth.n;
        for (unsigned j=i; j<_duration.n_elem; ++j)
        {
          _state.slice(j) = states;
          _loglik.at(j)   = -trans_mat.loglikelihood(states, n, _y.slice(j), _transition.slice(j));
        }
      }
      std::cout << "Deviance after smoothing round " << s << ": " << 2*arma::accu(_loglik) << std::endl;
      if (s > 0 && old_loglik - arma::accu(_loglik) < 0.0001) break;
      old_loglik = arma::accu(_loglik);
    }

    return arma::accu(_loglik);
  }

  double loglikelihood (void) const
  {
    return arma::accu(_loglik);
  }

  arma::cube parameters (void) const
  {
    return _parameters;
  }

  arma::cube state (void) const
  {
    return _state;
  }

  arma::cube transition (void) const
  {
    return _transition;
  }
};

struct particoal_lambda
{
  // struct particle above is really inefficient for updating parameters at a single window
  // and calculating likelihood. this does that much more efficiently by caching the transition matrix
  // and state variables.

  transition_matrix trans_mat;

  arma::vec lambda;

  arma::ucube _y, _mask;
  arma::umat _remap, _n;
  arma::vec  _duration, _loglik;
  arma::cube _parameters, _state, _transition;

  particoal_lambda (const unsigned num_pop, const unsigned num_ind)
  : trans_mat (num_pop, num_ind)
  {}

  arma::umat transitory_states (void)
  {
    // internal transitory states

    return trans_mat.t_states;
  }

  arma::umat coalescent_states (void)
  {
    // internal coalescent states

    return trans_mat.c_states;
  }

  arma::umat emission_classes (void)
  {
    // rows are starting states, columns are populations, values are haplotype counts
    // these classes correspond to columns in data "y" / elements of "z"
    // used in loglikelihood()

    return trans_mat.s_states;
  }

  arma::umat emission_states (void)
  {
    // rows are emission types (types of pairwise coalescences), columns are first/second of pair, values are population indices
    // these states correspond to rows in data "y" used in loglikelihood()

    return trans_mat.e_states;
  }

  arma::mat transition_probabilities (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_probabilities(remap, migr_mat, duration);
  }

  arma::mat transition_rates (const arma::uvec& remap, const arma::mat& migr_mat, const double& duration)
  {
    return trans_mat.transition_rates(remap, migr_mat, duration);
  }

  struct internal_nlopt_data
  {
    // these should really be const
    unsigned i;
    particoal_lambda& p;

    // these are mutable
    int        err;
    unsigned   num_par;
    arma::vec  start, lambda;
    arma::umat mapping;
    arma::mat  parameters, last_parameters, next_parameters, states;
    arma::uvec n;
    double     loglik;

    // settings
    bool verbose = true;
    bool smooth = false;

    static double filter (unsigned dim, const double* x, double* grad, void* f_data)
    {
      internal_nlopt_data* data = static_cast<internal_nlopt_data*>(f_data);

      unsigned i = data->i;

      data->n      = data->p._n.col(i);
      data->states = data->p._state.slice(i);

      if (arma::max(arma::vectorise(data->mapping)) > dim) Rcpp::stop("bad mask");

      data->parameters.zeros();
      for (unsigned j=0; j<data->mapping.n_elem; ++j)
        if (data->mapping.at(j) > 0)
          data->parameters.at(j) = exp(x[data->mapping.at(j)-1]);

      double lp = data->p.trans_mat.loglikelihood(data->states, data->n, data->p._y.slice(i), data->p._remap.col(i), data->parameters, data->p._duration.at(i));

      // random walk probabilities for pop sizes, etc
      // say we are right before a split. then mapping for next parameters might not defined.
      // if we are right after a split, we have 2 or more prior parameters
      //if (data->last_parameters.n_elem > 0)
      //{
      //  for (unsigned j=0; j<parameters.n_elem; ++j)
      //  {
      //    if (mapping.at(j) > 0)
      //    {
      //      lp += -lambda.at(j) * 
      //        std::pow(x[mapping.at(j)-1] - data->last_x[mapping.at(j)-1], 2);
      //    }
      //  }
      //}
      //if (data->next_parameters.n_elem > 0)
      //{
      //  for (unsigned j=0; j<parameters.n_elem; ++j)
      //  {
      //    if (mapping.at(j) > 0)
      //    {
      //      lp += -lambda.at(j) * 
      //        std::pow(x[mapping.at(j)-1] - data->next_x.at(j), 2);
      //    }
      //  }
      //}

      // if we're doing smoothing, we propagate the state vector forward and incorporate loglikelihood of future observations
      if (data->smooth)
        for (unsigned j=i+1; j<data->p._duration.n_elem; ++j)
          lp += data->p.trans_mat.loglikelihood(data->states, data->n, data->p._y.slice(j), data->p._transition.slice(j));

      return -lp;
    }

    internal_nlopt_data (particoal_lambda& p, const unsigned i, bool do_smooth, bool symmetric_migration) : i(i), p(p)
    {
      if (i >= p._duration.n_elem) Rcpp::stop("bad i");

      smooth = do_smooth;

      parameters = p._parameters.slice(i);

      // transition likelihood for parameters depends on subsequent and prior state
      //last_parameters = i > 0 ? 
      //  p._parameters.slice(i-1) : arma::zeros<arma::mat>(0,0);
      //next_parameters = i+1 < p._duration.n_elem ? 
      //  p._parameters.slice(i+1) : arma::zeros<arma::mat>(0,0);

      // find set of active parameters and linear indices thereof
      arma::uvec active = arma::unique(p._remap.col(i));
      mapping = arma::zeros<arma::umat>(p._remap.n_rows, p._remap.n_rows);
      num_par = 0;
      if (symmetric_migration)
      {
        for (auto j : active)
          for (auto k : active)
          {
            if (k >= j)
            {
              num_par += 1;
              mapping.at(k,j) = mapping.at(j,k) = num_par;
            }
          }
      } else {
        for (auto j : active)
          for (auto k : active)
          {
            num_par += 1;
            mapping.at(k,j) = num_par;
          }
      }

      // get starting values
      start = std::vector<double>(num_par);
      for (unsigned j=0; j<parameters.n_elem; ++j)
      {
        if (mapping.at(j) > 0)
        {
          start[mapping.at(j)-1] = log(parameters.at(j));
        }
      }

      nlopt_opt optimizer = nlopt_create(NLOPT_LN_BOBYQA, num_par);
      nlopt_set_min_objective(optimizer, filter, this);
      nlopt_set_xtol_rel(optimizer, 1e-4);
      err = nlopt_optimize(optimizer, &(start[0]), &loglik);

      if (err < 0)
      { // error; revert to initial state
        if (verbose) std::cout << "[nlopt] Failed with code: " << err << std::endl;
        parameters.print("Dumping parameters then reverting:");
        parameters = p._parameters.slice(i);
        for (unsigned j=0; j<parameters.n_elem; ++j)
          if (mapping.at(j) > 0)
            start[mapping.at(j)-1] = log(parameters.at(j));
      } else {
        for (unsigned j=0; j<parameters.n_elem; ++j)
          if (mapping.at(j) > 0)
            parameters.at(j) = exp(start[mapping.at(j)-1]);
      }

      // ensure that state/n are consistent with optimized parameters; and that loglik is just for the current observation
      smooth = false;
      loglik = filter(num_par, &(start[0]), nullptr, this);

      nlopt_destroy(optimizer);
    }
  };

  double estimate (arma::uvec n, const arma::ucube& y, const arma::umat& remap, const arma::cube& migr_mat, const arma::vec& duration, const unsigned num_smooth)
  {
    // calculate all transition matrices then loglikelihood
    // store n (starting lineages) for each window
    // also store _starting_ state for each window

    bool symmetric_migration = true;

    if(!(y.n_slices == duration.n_elem)) Rcpp::stop("invalid y");
    if(!(remap.n_cols == duration.n_elem)) Rcpp::stop("invalid remap");
    if(!(migr_mat.n_slices == duration.n_elem)) Rcpp::stop("invalid migr_mat");

    _y = y;
    _remap = remap;
    _parameters = migr_mat;
    _duration = duration;
    _n = arma::zeros<arma::umat>(n.n_elem, duration.n_elem);
    _loglik = arma::zeros<arma::vec>(duration.n_elem);
    _transition = arma::zeros<arma::cube>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans+trans_mat.num_coal, duration.n_elem);
    _state = arma::zeros<arma::cube>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_start, duration.n_elem);

    // filtering step
    arma::mat states = arma::eye<arma::mat>(trans_mat.num_trans+trans_mat.num_coal, trans_mat.num_trans);
    states = states.cols(trans_mat.s_map_t);
    for (unsigned i=0; i<_duration.n_elem; ++i)
    {
      // inputs to optimization
      _state.slice(i) = states;
      _n.col(i) = n;

      internal_nlopt_data filter (*this, i, false, symmetric_migration);

      // update running variables
      n = filter.n;
      states = filter.states;

      // outputs of optimization
      _loglik.at(i) = filter.loglik;
      _parameters.slice(i) = filter.parameters;
      _transition.slice(i) = trans_mat.transition_probabilities(_remap.col(i), _parameters.slice(i), _duration.at(i));
    }
    std::cout << "Deviance after filtering: " << 2.*arma::accu(_loglik) << std::endl;

    // smoothing steps
    double old_loglik;
    for (unsigned s=0; s<num_smooth; ++s)
    {
      for (unsigned i=_duration.n_elem; i>0; --i)
      {
        internal_nlopt_data smooth (*this, i-1, true, symmetric_migration);

        // outputs of optimization
        _loglik.at(i-1)        = smooth.loglik;
        _parameters.slice(i-1) = smooth.parameters;
        _transition.slice(i-1) = trans_mat.transition_probabilities(_remap.col(i-1), _parameters.slice(i-1), _duration.at(i-1));

        // propagate the new states forward and update loglikelihood of future timepoints
        states = smooth.states;
        n      = smooth.n;
        for (unsigned j=i; j<_duration.n_elem; ++j)
        {
          _state.slice(j) = states;
          _loglik.at(j)   = -trans_mat.loglikelihood(states, n, _y.slice(j), _transition.slice(j));
        }
      }
      std::cout << "Deviance after smoothing round " << s << ": " << 2*arma::accu(_loglik) << std::endl;
      if (s > 0 && old_loglik - arma::accu(_loglik) < 0.0001) break;
      old_loglik = arma::accu(_loglik);
    }

    return arma::accu(_loglik);
  }

  double loglikelihood (void) const
  {
    return arma::accu(_loglik);
  }

  arma::cube parameters (void) const
  {
    return _parameters;
  }

  arma::cube state (void) const
  {
    return _state;
  }

  arma::cube transition (void) const
  {
    return _transition;
  }
};

struct tally_tmrca
{
  unsigned num_pop, num_hap;
  arma::ucube e_counts;
  arma::uvec s_states, n_trees;
  arma::umat s_states_definition, e_states_definition;
  arma::uvec s_map_t;
  arma::umat c_map_e;

  bool collapse = true; //collapse clades that have equivalent population configuration
  
  //TODO factorize this
  tally_tmrca (const arma::cube& tmrca, const arma::umat& clades, const arma::uvec& populations, const arma::vec& breakpoints)
  {
    arma::uvec unique_populations = arma::unique(populations);

    num_pop = unique_populations.n_elem;
    num_hap = clades.n_rows;

    if (!(num_pop == unique_populations.max()+1 && clades.max() < populations.n_elem)) Rcpp::stop("invalid populations");
    if (!(breakpoints.at(0) == 0.0 && arma::all(arma::diff(breakpoints) > 0.0))) Rcpp::stop("invalid breakpoints");
    //if (!(clades.n_rows < clades.n_cols)) Rcpp::stop("invalid clades");

    transition_matrix trans_mat (num_pop, num_hap);
    e_states_definition = trans_mat.e_states;
    s_states_definition = trans_mat.s_states;
    s_map_t = trans_mat.s_map_t;
    c_map_e = trans_mat.c_map_e;
    n_trees = tmrca.n_slices * arma::ones<arma::uvec>(clades.n_cols);

    arma::mat windows (2, breakpoints.n_elem);
    for (unsigned i=0; i<breakpoints.n_elem-1; ++i)
    {
      windows.at(0,i) = breakpoints.at(i);
      windows.at(1,i) = breakpoints.at(i+1);
    }
    windows.at(0,windows.n_cols-1) = breakpoints.at(windows.n_cols-1);
    windows.at(1,windows.n_cols-1) = arma::datum::inf;

    e_counts = arma::zeros<arma::ucube>(trans_mat.e_states.n_rows, clades.n_cols, windows.n_cols);
    s_states = arma::zeros<arma::uvec>(clades.n_cols);

    // over pairs/trios/quads/whatever of haplotypes
    for (unsigned cl=0; cl<clades.n_cols; ++cl)
    {
      arma::uvec clade = clades.col(cl);

      // classify type of each possible MRCA (e.g. which populations coalesce)
      arma::uvec pair_to_e (clade.n_elem * (clade.n_elem-1) / 2); 
      unsigned k = 0;
      for (auto i : clade)
        for (auto j : clade)
          if (i < j)
          {
            unsigned p1 = populations.at(i), p2 = populations.at(j);
            arma::uvec e_state = 
              arma::find(trans_mat.e_states.col(0) == std::min(p1,p2) && 
                  trans_mat.e_states.col(1) == std::max(p1,p2), 1);
            if (e_state.n_elem == 0) Rcpp::stop("out of range in e_states");
            pair_to_e.at(k) = e_state[0];
            k++;
          }

      // classify clade into starting state (e.g. starting population configuration)
      arma::uvec s_state (trans_mat.s_states.n_cols, arma::fill::zeros);
      for (auto i : clade) s_state.at(populations.at(i))++;
      unsigned s = 0;
      bool match;
      for (; s<trans_mat.s_states.n_rows; ++s)
      {
        match = true;
        for (unsigned j=0; j<trans_mat.s_states.n_cols; ++j)
          match = match && trans_mat.s_states.at(s,j) == s_state.at(j);
        if (match) break;
      }
      if (!match) Rcpp::stop("out of range in s_states");
      s_states.at(cl) = s;

      // find first coalescent event and its emission state (e.g. between what populations)
      arma::uvec e_index (tmrca.n_slices);
      arma::vec c_time (tmrca.n_slices);
      for (unsigned i=0; i<tmrca.n_slices; ++i)
      {
        int min_i = -1;
        double min_t = arma::datum::inf;
        unsigned l = 0;
        for (auto j : clade)
        {
          for (auto k : clade)
          {
            if (j < k)
            {
              if (tmrca.at(j,k,i) < min_t) {min_t = tmrca.at(j,k,i); min_i = l;}
              l++;
            }
          }
        }
        if (min_i < 0 || min_i >= pair_to_e.n_elem) Rcpp::stop("out of range in tmrca");
        e_index.at(i) = pair_to_e.at(min_i);
        c_time.at(i) = min_t;
      }

      // classify coalescent events into time windows
      arma::uvec c_time_order = arma::stable_sort_index(c_time);
      unsigned w = 0;
      for (auto i : c_time_order)
      {
        if (!(c_time.at(i) >= windows.at(0,w) && c_time.at(i) < windows.at(1,w)))
        {
          do 
          { 
            w++; 
            if (c_time.at(i) >= windows.at(0,w) && c_time.at(i) < windows.at(1,w))
              break;
          } while (true);
        }
        e_counts.at(e_index.at(i), cl, w)++;
      }
    }

    // collapse from clades to starting configurations
    // for now, ignoring last window wherein every remaining lineage coalesces, although ultimately i could modify likelihood to use this info
    if (collapse)
    {
      arma::ucube e_counts_collapse (trans_mat.e_states.n_rows, trans_mat.s_states.n_rows, windows.n_cols-1, arma::fill::zeros);
      for (unsigned i=0; i<trans_mat.e_states.n_rows; ++i)
        for (unsigned j=0; j<clades.n_cols; ++j)
          for (unsigned k=0; k<windows.n_cols-1; ++k)
            e_counts_collapse.at(i,s_states.at(j),k) += e_counts.at(i,j,k);
      e_counts = e_counts_collapse;

      arma::uvec n_trees_collapse (trans_mat.s_states.n_rows, arma::fill::zeros);
      for (unsigned i=0; i<clades.n_cols; ++i)
        n_trees_collapse.at(s_states.at(i)) += n_trees.at(i);
      n_trees = n_trees_collapse;
    }
  }

  arma::ucube get_e_counts (void) const { return e_counts; };
  arma::uvec get_s_states (void) const { return s_states; };
  arma::uvec get_n_trees (void) const { return n_trees; };
  arma::umat get_e_states_definition (void) const { return e_states_definition; };
  arma::umat get_s_states_definition (void) const { return s_states_definition; };
  arma::uvec get_s_map_t (void) const { return s_map_t; };
  arma::umat get_c_map_e (void) const { return c_map_e; };
};

RCPP_EXPOSED_CLASS_NODECL(particle)

RCPP_MODULE(particle) {
  using namespace Rcpp;
  class_<decoder>("decoder")
    .constructor<unsigned, unsigned, bool>()
    .method("emission_classes", &decoder::emission_classes)
    .method("emission_states", &decoder::emission_states)
    .method("transitory_states", &decoder::transitory_states)
    .method("coalescent_states", &decoder::coalescent_states)
    .method("map_coalescent_states_to_emission_states", &decoder::map_coalescent_states_to_emission_states)
    .method("map_emission_classes_to_transitory_states", &decoder::map_emission_classes_to_transitory_states)
    .method("symbolic_transition_rates", &decoder::symbolic_transition_rates)
    .method("transition_probabilities", &decoder::transition_probabilities)
    .method("transition_rates", &decoder::transition_rates)
    .method("migration_function", &decoder::migration_function)
    .method("migration_operator", &decoder::migration_operator)
    .method("migration_function_penalty", &decoder::migration_function_penalty)
    .method("debug", &decoder::debug)
    .method("deviance", &decoder::deviance)
    .method("simulate", &decoder::simulate)
    .method("initial_states", &decoder::initial_states)
    ;
  class_<particle>("particle")
    .constructor<unsigned, unsigned>()
    .method("emission_classes", &particle::emission_classes)
    .method("emission_states", &particle::emission_states)
    .method("transitory_states", &particle::transitory_states)
    .method("coalescent_states", &particle::coalescent_states)
    .method("map_coalescent_states_to_emission_states", &particle::map_coalescent_states_to_emission_states)
    .method("map_emission_classes_to_transitory_states", &particle::map_emission_classes_to_transitory_states)
    .method("transition_probabilities", &particle::transition_probabilities)
    .method("symbolic_transition_rates", &particle::symbolic_transition_rates)
    .method("transition_rates", &particle::transition_rates)
    .method("loglikelihood", &particle::loglikelihood)
    .method("simulate", &particle::simulate)
    ;
  class_<particulate>("particulate")
    .constructor<unsigned, unsigned>()
    .method("emission_classes", &particulate::emission_classes)
    .method("emission_states", &particulate::emission_states)
    .method("transitory_states", &particulate::transitory_states)
    .method("coalescent_states", &particulate::coalescent_states)
    .method("transition_probabilities", &particulate::transition_probabilities)
    .method("transition_rates", &particulate::transition_rates)
    .method("loglikelihood", &particulate::loglikelihood)
    .method("update", &particulate::update)
    .method("init", &particulate::init)
    .method("parameters", &particulate::parameters)
    .method("state", &particulate::state)
    .method("transition", &particulate::transition)
    ;
  class_<particoal>("particoal")
    .constructor<unsigned, unsigned>()
    .method("emission_classes", &particoal::emission_classes)
    .method("emission_states", &particoal::emission_states)
    .method("transitory_states", &particoal::transitory_states)
    .method("coalescent_states", &particoal::coalescent_states)
    .method("transition_probabilities", &particoal::transition_probabilities)
    .method("transition_rates", &particoal::transition_rates)
    .method("loglikelihood", &particoal::loglikelihood)
    .method("estimate", &particoal::estimate)
    .method("parameters", &particoal::parameters)
    .method("state", &particoal::state)
    .method("transition", &particoal::transition)
    ;
  class_<particoal_lambda>("particoal_lambda")
    .constructor<unsigned, unsigned>()
    .method("emission_classes", &particoal_lambda::emission_classes)
    .method("emission_states", &particoal_lambda::emission_states)
    .method("transitory_states", &particoal_lambda::transitory_states)
    .method("coalescent_states", &particoal_lambda::coalescent_states)
    .method("transition_probabilities", &particoal_lambda::transition_probabilities)
    .method("transition_rates", &particoal_lambda::transition_rates)
    .method("loglikelihood", &particoal_lambda::loglikelihood)
    .method("estimate", &particoal_lambda::estimate)
    .method("parameters", &particoal_lambda::parameters)
    .method("state", &particoal_lambda::state)
    .method("transition", &particoal_lambda::transition)
    ;
  class_<tally_tmrca>("tally_tmrca")
    .constructor<arma::cube, arma::umat, arma::uvec, arma::vec>()
    .method("get_e_counts", &tally_tmrca::get_e_counts)
    .method("get_s_states", &tally_tmrca::get_s_states)
    .method("get_n_trees", &tally_tmrca::get_n_trees)
    .method("get_s_states_definition", &tally_tmrca::get_s_states_definition)
    .method("get_e_states_definition", &tally_tmrca::get_e_states_definition)
    .method("get_s_map_t", &tally_tmrca::get_s_map_t)
    .method("get_c_map_e", &tally_tmrca::get_c_map_e)
    ;
}

