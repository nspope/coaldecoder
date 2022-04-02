#ifndef _TRANSITION_RATE_MATRIX
#define _TRANSITION_RATE_MATRIX

#include <RcppArmadillo.h> 
#include <vector>
#include <string>
#include <algorithm>

struct TrioTransitionRates
{
  /*
   *  Transition rates between trios under an approximate structured
   *  coalescent, wherein only one transition can happen instantaneously. This
   *  approximation is very accurate if rates are small, and makes the rate
   *  matrix sparse.
   *
   *  Transition rates are parameterized by demographic parameters M, where
   *  M[i,j] is the backwards-in-time migration of lineages from population i
   *  to population j; and M[i,i] is the haploid effective population size of
   *  population i.
   */

  const std::string prefix = "[TrioTransitionRates] ";
  const bool check_valid = true;

  const unsigned P; //number of populations
  const arma::uvec::fixed<3> S; //number of 3,2,1-lineage states
  const arma::mat M; //demographic parameters
  const arma::sp_mat X; //rate matrix

  TrioTransitionRates (const arma::mat& _M, const bool _check_valid = true)
    : check_valid (_check_valid)
    , P (_M.n_rows)
    , S ({power(P, 3), 3*power(P, 2), 3*P})
    , M (_M)
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

  std::string initial_state (arma::uvec::fixed<3> initial) const
  {
    /*
     *  Nonredundant population labelling of trio, "{u_1, u_2, u_3}"
     */

    std::sort(initial.begin(), initial.end());
    std::string out = "{";
    out += std::to_string(initial[0]) + ",";
    out += std::to_string(initial[1]) + ",";
    out += std::to_string(initial[2]) + "}";
    return out;
  }

  std::vector<std::string> initial_states (void) const
  {
    std::vector<std::string> out;
    arma::uvec::fixed<3> u;
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=u[1]; u[2]<P; ++u[2])
        {
          out.emplace_back(initial_state(u));
        }
      }
    }

    // assert that string array is sorted
    if (!std::is_sorted(out.begin(), out.end()))
    {
      Rcpp::stop(prefix + " Initial states incorrectly named");
    }

    return out;
  }

  std::string emission_state (const arma::uvec::fixed<3>& initial, const arma::uvec::fixed<3>& absorbing) const
  {
    /*
     *  2-lineage states:
     *  [P, ., .] ==> ((u[0], u[1])x, u[2])
     *  [., P, .] ==> ((u[1], u[2])x, u[0])
     *  [., ., P] ==> ((u[2], u[0])x, u[1])
     * 
     *  3-lineage states:
     *  [P, P, .] ==> ((u[0], u[1]), u[2])x
     *  [., P, P] ==> ((u[1], u[2]), u[0])x
     *  [P, ., P] ==> ((u[2], u[0]), u[1])x
     * 
     *  Where 'u' is the starting state. The ordering of 'u' within the inner
     *  clade is always (min, max).
     */

    unsigned surviving_lineages = arma::accu(absorbing != P);
    if (arma::any(absorbing > P) || surviving_lineages == 0 || surviving_lineages == 3)
    {
      Rcpp::stop(prefix + "Absorbing state out-of-range");
    }

    if (arma::any(initial >= P))
    {
      Rcpp::stop(prefix + "Initial state out-of-range");
    }

    arma::uvec::fixed<3> u;
    if (surviving_lineages == 1)
    {
      arma::uvec v = arma::find(absorbing != P);
      u[0] = initial[(v[0]+1) % 3];
      u[1] = initial[(v[0]+2) % 3];
      u[2] = initial[(v[0]+0) % 3];
    } else {
      arma::uvec v = arma::find(absorbing == P);
      u[0] = initial[(v[0]+0) % 3];
      u[1] = initial[(v[0]+1) % 3];
      u[2] = initial[(v[0]+2) % 3];
    }
    u.head(2) = arma::sort(u.head(2));
    std::string out = 
      surviving_lineages == 2 ? "t1::((" : "t2::((";
    out += std::to_string(u[0]) + ",";
    out += std::to_string(u[1]) + "),";
    out += std::to_string(u[2]) + ")";
    return out;
  }

  std::vector<std::string> emission_states (void) const
  {
    std::vector<std::string> out;
    arma::uvec::fixed<3> u;

    // 2-lineage states: "((u[0],u[1]),u[2]).t1" 
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=0; u[2]<P; ++u[2])
        {
          arma::uvec::fixed<3> v = {P, 0, 0}; 
          out.emplace_back(emission_state(u, v));
        }
      }
    }

    // 1-lineage states: "((u[0],u[1]),u[2]).t2"
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=0; u[2]<P; ++u[2])
        {
          arma::uvec::fixed<3> v = {P, P, 0}; 
          out.emplace_back(emission_state(u, v));
        }
      }
    }

    // assert that string array is sorted
    if (!std::is_sorted(out.begin(), out.end()))
    {
      Rcpp::stop(prefix + " Emission states incorrectly named");
    }

    return out;
  }

  arma::umat initial_to_states (void) const
  {
    /*
     *  Map from initial population labelling onto first compatible state
     */

    std::vector<std::string> initial_names = initial_states();

    std::vector<arma::uword> col_indices;
    std::vector<arma::uword> row_indices;
    col_indices.reserve(R::choose(P + 2, P - 2));
    row_indices.reserve(R::choose(P + 2, P - 2));

    arma::uvec::fixed<3> u;
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=u[1]; u[2]<P; ++u[2])
        {
          auto initial_index = std::lower_bound(
            initial_names.begin(), 
            initial_names.end(), 
            initial_state(u)
          );
          col_indices.emplace_back(initial_index - initial_names.begin());
          row_indices.emplace_back(linear_index(u));
        }
      }
    }

    arma::umat out = arma::join_vert(
      arma::urowvec(col_indices),
      arma::urowvec(row_indices)
    );
    return out;
  }

  arma::umat states_to_emission (void) const
  {
    /*
     *  Map from coalescent states onto emission states
     */

    std::vector<std::string> initial_names = initial_states();
    std::vector<std::string> emission_names = emission_states();

    std::vector<arma::uword> col_indices;
    std::vector<arma::uword> row_indices;
    std::vector<arma::uword> out_indices;
    col_indices.reserve(S[1] + S[2]);
    row_indices.reserve(S[1] + S[2]);
    out_indices.reserve(S[1] + S[2]);

    arma::uvec::fixed<3> u; // non-redundant starting configuration
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=u[1]; u[2]<P; ++u[2])
        {
          auto initial_index = std::lower_bound(
            initial_names.begin(), 
            initial_names.end(), 
            initial_state(u)
          );
          if (initial_index == initial_names.end())
          {
            Rcpp::stop(prefix + " Initial state out of range");
          }
          arma::uvec::fixed<3> v;

          // coalescent states with two surviving lineages
          for (unsigned i=0; i<3; ++i)
          {
            v[i] = P;
            unsigned j = (i + 1) % 3;
            unsigned k = (i + 2) % 3;
            for (v[j] = 0; v[j] < P; ++v[j])
            {
              for (v[k] = 0; v[k] < P; ++v[k])
              {
                auto emission_index = std::lower_bound(
                  emission_names.begin(), 
                  emission_names.end(), 
                  emission_state(u, v)
                );
                if (emission_index == emission_names.end())
                {
                  std::cout << "OOR: " << emission_state(u, v) << std::endl;
                  Rcpp::stop(prefix + " Emission state out-of-range");
                }
                col_indices.emplace_back(initial_index - initial_names.begin());
                row_indices.emplace_back(linear_index(v));
                out_indices.emplace_back(emission_index - emission_names.begin());
              }
            }
          }

          // coalescent states with one surviving lineage
          for (unsigned i=0; i<3; ++i)
          {
            unsigned j = (i + 1) % 3;
            unsigned k = (i + 2) % 3;
            v[j] = P;
            v[k] = P;
            for (v[i]=0; v[i]<P; ++v[i])
            {
              auto emission_index = std::lower_bound(
                emission_names.begin(), 
                emission_names.end(), 
                emission_state(u, v)
              );
              if (emission_index == emission_names.end())
              {
                std::cout << "OOR: " << emission_state(u, v) << std::endl;
                Rcpp::stop(prefix + " Emission state out-of-range");
              }
              col_indices.emplace_back(initial_index - initial_names.begin());
              row_indices.emplace_back(linear_index(v));
              out_indices.emplace_back(emission_index - emission_names.begin());
            }
          }
        }
      }
    }

    arma::umat out = arma::join_vert(
        arma::urowvec(col_indices),
        arma::urowvec(row_indices),
        arma::urowvec(out_indices)
    );

    return out;
  }

  arma::umat emission_to_initial (void) const
  {
    /*
     *  Indices of initial state associated with emission, and which coalescence event it is
     */

    std::vector<std::string> initial_names = initial_states();
    std::vector<std::string> emission_names = emission_states();

    arma::umat out (2, emission_names.size());
    arma::uvec::fixed<3> u;

    // 2-lineage states
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=0; u[2]<P; ++u[2])
        {
          arma::uvec::fixed<3> v = {P, 0, 0}; 
          auto initial_index = std::lower_bound(
              initial_names.begin(), 
              initial_names.end(), 
              initial_state(u)
          );
          auto emission_index = std::lower_bound(
              emission_names.begin(), 
              emission_names.end(), 
              emission_state(u, v)
          );
          out.at(0, emission_index - emission_names.begin()) = 
            initial_index - initial_names.begin();
          out.at(1, emission_index - emission_names.begin()) = 2;
        }
      }
    }

    // 1-lineage states: "((u[0],u[1]),u[2]).t2"
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=u[0]; u[1]<P; ++u[1])
      {
        for (u[2]=0; u[2]<P; ++u[2])
        {
          arma::uvec::fixed<3> v = {P, P, 0}; 
          auto initial_index = std::lower_bound(
              initial_names.begin(), 
              initial_names.end(), 
              initial_state(u)
          );
          auto emission_index = std::lower_bound(
              emission_names.begin(), 
              emission_names.end(), 
              emission_state(u, v)
          );
          out.at(0, emission_index - emission_names.begin()) = 
            initial_index - initial_names.begin();
          out.at(1, emission_index - emission_names.begin()) = 1;
        }
      }
    }

    return out;
  }

  arma::sp_mat build_matrix (void) const
  {
    /*
     *  Return rate matrix
     */

    if (M.n_rows != P || M.n_cols != P)
    {
      Rcpp::stop(prefix + "Parameter matrix dimensions do not equal populations");
    }

    if (check_valid)
    {
      if (arma::any(M.diag() == 0.0) || arma::any(arma::vectorise(M) < 0.0))
      {
        Rcpp::stop(prefix + "Invalid parameter matrix");
      }
    }

    arma::uword storage_bound =  2 * (
      S[0] * (P-1) * 3 + // (P-1) * 3 migrations per 3-lineage
      S[1] + // first coalescent transitions
      S[1] * (P-1) * 2 + // (P-1) * 2 migrations per 2-lineage
      S[2] + // second coalescent transitions
      S[2] * (P-1) // (P-1) * 1 migrations per 1-lineage
     );

    std::vector<arma::uword> rows, cols;
    std::vector<double> values;
    rows.reserve(storage_bound);
    cols.reserve(storage_bound);
    values.reserve(storage_bound);

    // three lineages
    arma::uvec::fixed<3> u;
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=0; u[1]<P; ++u[1])
      {
        for (u[2]=0; u[2]<P; ++u[2])
        {
          arma::uword source_index = linear_index(u);

          // lineages migrate
          for (unsigned i=0; i<3; ++i)
          {
            arma::uvec::fixed<3> _u = u;
            for (_u[i]=0; _u[i]<P; ++_u[i])
            {
              if (u[i] != _u[i])
              {
                arma::uword dest_index = linear_index(_u);
                push_back(rows, cols, values, source_index, dest_index, M.at(u[i],_u[i]));
                push_back(rows, cols, values, source_index, source_index, -M.at(u[i],_u[i]));
                //rows.push_back(source_index);
                //cols.push_back(dest_index);
                //values.push_back(M(u[i],_u[i]));
                //rows.push_back(source_index);
                //cols.push_back(source_index);
                //values.push_back(-1.0 * M(u[i],_u[i]));
              }
            }
          }

          // lineages coalesce
          for (unsigned i=0; i<3; ++i)
          {
            unsigned j = (i + 1) % 3;
            if (u[i] == u[j])
            {
              arma::uvec _u = u; _u[i] = P;
              arma::uword dest_index = linear_index(_u);
              push_back(rows, cols, values, source_index, dest_index, 1.0/M.at(u[i],u[j]));
              push_back(rows, cols, values, source_index, source_index, -1.0/M.at(u[i],u[j]));
              //rows.push_back(source_index);
              //cols.push_back(dest_index);
              //values.push_back(1.0 / M(u[i],u[j]));
              //rows.push_back(source_index);
              //cols.push_back(source_index);
              //values.push_back(-1.0 / M(u[i],u[j]));
            }
          }
        }
      }
    }

    // two lineages
    for (unsigned i=0; i<3; ++i)
    {
      u[i] = P;
      unsigned j = (i + 1) % 3;
      unsigned k = (i + 2) % 3;

      for (u[j] = 0; u[j] < P; ++u[j])
      {
        for (u[k] = 0; u[k] < P; ++u[k])
        {
          arma::uword source_index = linear_index(u);

          // lineages migrate
          arma::uvec::fixed<2> v = {j, k};
          for (auto l : v)
          {
            arma::uvec::fixed<3> _u = u;
            for (_u[l] = 0; _u[l] < P; ++_u[l])
            {
              if (_u[l] != u[l])
              {
                arma::uword dest_index = linear_index(_u);
                push_back(rows, cols, values, source_index, dest_index, M.at(u[l],_u[l]));
                push_back(rows, cols, values, source_index, source_index, -M.at(u[l],_u[l]));
                //rows.push_back(source_index);
                //cols.push_back(dest_index);
                //values.push_back(M(u[l],_u[l]));
                //rows.push_back(source_index);
                //cols.push_back(source_index);
                //values.push_back(-1.0 * M(u[l],_u[l]));
              }
            }
          }

          // lineages coalesce
          if (u[j] == u[k])
          {
            arma::uvec::fixed<3> _u = u; _u[j] = P;
            arma::uword dest_index = linear_index(_u);
            push_back(rows, cols, values, source_index, dest_index, 1.0/M.at(u[j],u[k]));
            push_back(rows, cols, values, source_index, source_index, -1.0/M.at(u[j],u[k]));
            //rows.push_back(source_index);
            //cols.push_back(dest_index);
            //values.push_back(1.0 / M(u[j],u[k]));
            //rows.push_back(source_index);
            //cols.push_back(source_index);
            //values.push_back(-1.0 / M(u[j],u[k]));
          }
        }
      }
    }

    // one lineage
    for (unsigned i=0; i<3; ++i)
    {
      unsigned j = (i + 1) % 3;
      unsigned k = (i + 2) % 3;
      u[j] = P;
      u[k] = P;

      for (u[i]=0; u[i]<P; ++u[i])
      {
        arma::uword source_index = linear_index(u);
        arma::uvec::fixed<3> _u = u;
        for (_u[i]=0; _u[i]<P; ++_u[i])
        {
          if (u[i] != _u[i])
          {
            arma::uword dest_index = linear_index(_u);
            push_back(rows, cols, values, source_index, dest_index, M.at(u[i],_u[i]));
            push_back(rows, cols, values, source_index, source_index, -M.at(u[i],_u[i]));
            //rows.push_back(source_index);
            //cols.push_back(dest_index);
            //values.push_back(M(u[i],_u[i]));
            //rows.push_back(source_index);
            //cols.push_back(source_index);
            //values.push_back(-1.0 * M(u[i],_u[i]));
          }
        }
      }
    }

    /*
     *  3-lineage states: (P * P * P)
     *     [a,b,c] where a,b,c < P 
     *  2-lineage states: (3 * P * P)
     *     [P,a,b] where a,b < P ==> (0 x 1) have coalesced
     *     [a,P,b] where a,b < P ==> (1 x 2) have coalesced
     *     [a,b,P] where a,b < P ==> (2 x 0) have coalesced
     *  1-lineage states: (3 * P)
     *     [P,P,a] where a < P ==> (0 x 1) x 2
     *     [a,P,P] where a < P ==> (1 x 2) x 0
     *     [P,a,P] where a < P ==> (2 x 0) x 1
     *  [P,P,P] is unvisited
     */

    if (values.size() != storage_bound) 
    {
      Rcpp::stop(prefix + " Storage bound violated");
    }

    arma::umat locations = arma::join_vert(
        arma::urowvec(rows),
        arma::urowvec(cols)
    );
    arma::sp_mat out (true, locations, arma::vec(values), arma::accu(S), arma::accu(S));
    return out;
  }

  arma::mat reverse_differentiate (const arma::sp_mat& gradient) const
  {
    /*
     *  Use chain rule in reverse to map gradient back onto parameter matrix
     */

    if (gradient.n_rows != X.n_rows || gradient.n_cols != X.n_cols)
    {
      Rcpp::stop(prefix + "Gradient dimensions do not equal matrix dimensions");
    }

    arma::mat out (P, P, arma::fill::zeros);

    // three lineages
    arma::uvec::fixed<3> u;
    for (u[0]=0; u[0]<P; ++u[0])
    {
      for (u[1]=0; u[1]<P; ++u[1])
      {
        for (u[2]=0; u[2]<P; ++u[2])
        {
          arma::uword source_index = linear_index(u);

          // lineages migrate
          for (unsigned i=0; i<3; ++i)
          {
            arma::uvec::fixed<3> _u = u;
            for (_u[i]=0; _u[i]<P; ++_u[i])
            {
              if (u[i] != _u[i])
              {
                arma::uword dest_index = linear_index(_u);
                out.at(u[i],_u[i]) += gradient(source_index, dest_index);
                out.at(u[i],_u[i]) += -gradient(source_index, source_index);
              }
            }
          }

          // lineages coalesce
          for (unsigned i=0; i<3; ++i)
          {
            unsigned j = (i + 1) % 3;
            if (u[i] == u[j])
            {
              arma::uvec _u = u; _u[i] = P;
              arma::uword dest_index = linear_index(_u);
              out.at(u[i],u[j]) += -1.0/std::pow(M.at(u[i],u[j]), 2) * gradient(source_index, dest_index);
              out.at(u[i],u[j]) += 1.0/std::pow(M.at(u[i],u[j]), 2) * gradient(source_index, source_index);
            }
          }
        }
      }
    }

    // two lineages
    for (unsigned i=0; i<3; ++i)
    {
      u[i] = P;
      unsigned j = (i + 1) % 3;
      unsigned k = (i + 2) % 3;

      for (u[j] = 0; u[j] < P; ++u[j])
      {
        for (u[k] = 0; u[k] < P; ++u[k])
        {
          arma::uword source_index = linear_index(u);

          // lineages migrate
          arma::uvec::fixed<2> v = {j, k};
          for (auto l : v)
          {
            arma::uvec::fixed<3> _u = u;
            for (_u[l] = 0; _u[l] < P; ++_u[l])
            {
              if (_u[l] != u[l])
              {
                arma::uword dest_index = linear_index(_u);
                out.at(u[l],_u[l]) += gradient(source_index, dest_index);
                out.at(u[l],_u[l]) += -gradient(source_index, source_index);
              }
            }
          }

          // lineages coalesce
          if (u[j] == u[k])
          {
            arma::uvec::fixed<3> _u = u; _u[j] = P;
            arma::uword dest_index = linear_index(_u);
            out.at(u[j],u[k]) += -1.0/std::pow(M.at(u[j],u[k]), 2) * gradient(source_index, dest_index);
            out.at(u[j],u[k]) += 1.0/std::pow(M.at(u[j],u[k]), 2) * gradient(source_index, source_index);
          }
        }
      }
    }

    // one lineage
    for (unsigned i=0; i<3; ++i)
    {
      unsigned j = (i + 1) % 3;
      unsigned k = (i + 2) % 3;
      u[j] = P;
      u[k] = P;

      for (u[i]=0; u[i]<P; ++u[i])
      {
        arma::uword source_index = linear_index(u);
        arma::uvec::fixed<3> _u = u;
        for (_u[i]=0; _u[i]<P; ++_u[i])
        {
          if (u[i] != _u[i])
          {
            arma::uword dest_index = linear_index(_u);
            out.at(u[i],_u[i]) += gradient(source_index, dest_index);
            out.at(u[i],_u[i]) += -gradient(source_index, source_index);
          }
        }
      }
    }

    return out;
  }
};

#endif
