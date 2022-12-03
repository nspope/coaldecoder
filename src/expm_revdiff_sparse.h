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
#ifndef _EXPM_REVDIF_SPARSE
#define _EXPM_REVDIF_SPARSE

#include <RcppArmadillo.h> 
#include <map>
#include <string>
#include <random>

// C++ implementation based off of:
// https://github.com/scipy/scipy/blob/v1.8.0/scipy/sparse/linalg/_expm_multiply.py#L56-L142
// <scipy.linalg.onenormest> TODO
// Matrix::onenormest TODO

using namespace arma;

struct OneNormEst
{
  /*
   *  Calculate 1-norm of power of sparse matrix, exactly or approximately
   *  TODO: approximation seems to cause problems occasionally. check against scipy implementation
   */

  const std::string prefix = "[OneNormEst] ";
  const double eps = 4 * datum::eps;

  mutable std::mt19937 _rng;
  mutable std::uniform_real_distribution<> _randu;

  vec _v;
  vec _w;
  uword _iter;
  double _est;

  OneNormEst (const sp_mat& A, const int p=0, const bool approximate=true) :
    _rng(1024), _randu(0.0, 1.0)
  {
    //compute exact if !approximate or condition met?
    _est = approximate ? approximate_norm(A, p) : exact_norm(A, p);
  }

  double exact_norm (const sp_mat& A, const int p)
  {
    mat I = eye(A.n_rows, A.n_rows);
    for (unsigned i = 0; i<p; ++i)
    {
      I *= A;
    }
    return norm(I, 1);
  }

  mat A_x (const sp_mat& A, mat X, const int p)
  {
    for (unsigned i=0; i<p; ++i)
    {
      X = A * X;
    }
    return X;
  }

  mat At_x (const sp_mat& A, mat X, const int p)
  {
    for (unsigned i=0; i<p; ++i)
    {
      X = arma::trans(A) * X;
    }
    return X;
  }

  // ------- equivalent to `onenormest` in python::scipy ------- //
  // 07 JUN 22 -- I think this is working, gives equivalent results to scipy.
  //   However there are no resamples beyond the first. Should test this.

  double approximate_norm (const sp_mat& A, const int p, const unsigned t = 2, const unsigned iter_max = 5)
  {
    // Higham and Tisseur 2000, Algorithm 2.4
    
    uword n = A.n_rows;

    if (t >= n)
    {
      Rcpp::stop(prefix + "t >= n, use exact norm instead");
    }

    if (A.n_cols != n)
    {
      Rcpp::stop(prefix + "Matrix must be square");
    }

    if (iter_max < 2)
    {
      Rcpp::stop(prefix + "At least two iterations needed");
    }

    if (t < 1)
    {
      Rcpp::stop(prefix + "At least one column needed");
    }

    unsigned nmults = 0;
    unsigned nresamples = 0;
    // "We now explain our choice of starting matrix. We take the first
    // column of X to be the vector of 1s [...] This has the advantage that
    // for a matrix with nonnegative elements the algorithm converges
    // with an exact estimate on the second iteration, and such matrices
    // arise in applications [...]"
    mat X = ones(n, t); 
    mat Y = zeros(n, 0);
    mat Z = zeros(n, 0);
    mat S = zeros(n, t);
    mat S_old = zeros(n, 0);
    vec h = zeros(0);
    vec v;
    vec w;

    double est, est_old;

    // "The remaining columns are chosen as rand{-1,1},
    // with a check for and correction of parallel columns,
    // exactly as for S in the body of the algorithm."
    if (t > 1)
    {
      for (uword i=1; i<t; ++i)
      {
        // These are technically initial samples, not resamples,
        // so the resampling count is not incremented.
        _resample_column(i, X);
      }
      for (uword i=0; i<t; ++i)
      {
        do
        {
          if (_column_needs_resampling(i, X, Y))
          {
            _resample_column(i, X);
            nresamples++;
          } else {
            break;
          }
        } while (true);
      }
    }
    // "Choose starting matrix X with columns of unit 1-norm."
    X /= double(n);

    est_old = 0.0;
    int k = 1;
    uword ind_best;
    uvec ind = zeros<uvec>(0);
    // "indices of used unit vectors e_j"
    uvec ind_hist = zeros<uvec>(0);

    do
    {
      Y = A_x(A, X, p);
      nmults++;
      rowvec mags = _sum_abs_axis0(Y);
      est = mags.max();
      uword best_j = mags.index_max();
      if (est > est_old || k == 2)
      {
        if (k >= 2)
        {
          ind_best = ind(best_j);
        }
        w = Y.col(best_j);
      }
      // (1)
      if (k >= 2 && est <= est_old)
      {
        est = est_old;
        break;
      }
      est_old = est;
      S_old = S;
      if (k > iter_max) break;
      S = _sign_round_up(Y);
      Y.set_size(n, 0);
      // (2)
      if (_every_col_of_X_is_parallel_to_a_col_of_Y(S, S_old)) break;
      if (t > 1)
      {
        // "Ensure that no column of S is parallel to another column of S
        // or to a column of S_old by replacing columns of S by rand{-1,1}."
        for (uword i=0; i<t; ++i)
        {
          do
          {
            if (_column_needs_resampling(i, S, S_old))
            {
              _resample_column(i, S);
              nresamples++;
            } else {
              break;
            }
          } while (true);
        }
      }
      S_old.set_size(n, 0);
      // (3)
      Z = At_x(A, S, p);
      nmults++;
      h = _max_abs_axis1(Z);
      Z.set_size(n, 0);
      // (4)
      if (k >= 2 && h.max() == h(ind_best)) break;
      // "Sort h so that h_first >= ... >= h_last
      // and re-order ind correspondingly."
      //
      // Later on, we will need at most t+len(ind_hist) largest
      // entries, so drop the rest
      ind = stable_sort_index(h, "descend");
      ind = ind.head(t + ind_hist.n_elem);
      h.set_size(0);
      if (t > 1)
      {
        // (5)
        // Break if the most promising t vectors have been visited already.
        if (_set_difference(ind.head(t), ind_hist).n_elem == 0) break;
        // Put the most promising unvisited vectors at the front of the list
        // and put the visited vectors at the end of the list.
        // Preserve the order of the indices induces by the ordering of h.
        ind = join_vert(_set_difference(ind, ind_hist), _set_intersection(ind, ind_hist));
      }
      for (uword j=0; j<t; ++j)
      {
        X.col(j) = _elementary_vector(n, ind(j));
      }

      uvec new_ind = _set_difference(ind.head(t), ind_hist);
      ind_hist = join_vert(ind_hist, new_ind);
      k++;
    } while (true);

    v = _elementary_vector(n, ind_best); 

    _iter = k;
    _v = v;
    _w = w;

    return est;
  }

  // ------ internal ------ //

  vec _elementary_vector (const uword n, const uword i) const
  {
    vec out = zeros(n);
    out.at(i) = 1;
    return out;
  }

  bool _column_needs_resampling (const uword i, const mat& X, const mat& Y) const
  {
    uword n = X.n_rows;
    uword t = X.n_cols;
    bool out = false;
    for (uword j=0; j<i; ++j)
    {
      out = out || _vectors_are_parallel(X.col(i), X.col(j));
    }
    if (out) return out;
    for (uword j=0; j<Y.n_cols; ++j)
    {
      out = out || _vectors_are_parallel(X.col(i), Y.col(j));
    }
    return out;
  }

  vec _random_uniform (const uword n) const
  {
    vec out (n);
    for (uword i=0; i<n; ++i)
    {
      out[i] = _randu(_rng);
    }
    return out;
  }

  bool _vectors_are_parallel (const vec& v, const vec& w) const
  {
    return dot(v, w) == double(v.n_elem);
  }

  void _resample_column (const uword i, mat& X) const
  {
    //X.col(i) = sign(randu(X.n_rows) - 0.5);
    X.col(i) = sign(_random_uniform(X.n_rows) - 0.5);
  }

  mat _sign_round_up (arma::mat X) const
  {
    X.replace(0.0, 1.0);
    X /= abs(X);
    return X;
  }

  bool _every_col_of_X_is_parallel_to_a_col_of_Y (const mat& X, const mat& Y) const
  {
    for (uword i=0; i<X.n_cols; ++i)
    {
      for (uword j=0; j<Y.n_cols; ++j)
      {
        if (!_vectors_are_parallel(X.col(i), Y.col(j))) return false;
      }
    }
    return true;
  }

  rowvec _sum_abs_axis0 (const arma::mat& X)
  {
    return sum(abs(X), 0);
  }

  vec _max_abs_axis1 (const arma::mat& X)
  {
    vec out (X.n_rows);
    for (uword i=0; i<X.n_rows; ++i)
    {
      out.at(i) = abs(X.row(i)).max();
    }
    return out;
  }

  uvec _set_intersection (const uvec& A, const uvec& B)
  {
    // return in order that elements occur in A
    uvec Ai, Bi, C;
    intersect(C, Ai, Bi, A, B);
    Ai = sort(Ai);
    return A.elem(Ai);
  }

  uvec _set_difference (const uvec& A, const uvec& B)
  {
    // return in order that elements occur in A
    uvec Ai, Bi, C;
    intersect(C, Ai, Bi, A, B);
    uvec Ad = zeros<uvec>(A.n_elem);
    Ad.elem(Ai) += 1;
    return A.elem(find(Ad == 0));
  }
};

struct SparseMatrixExponentialMultiply
{
  const std::string prefix = "[SparseMatrixExponentialMultiply] ";
  const bool verbose = false;

  //private:
  std::map<int, double> _theta = {
    {1,2.29e-16},{2,2.58e-8},{3,1.39e-5},{4,3.40e-4},{5,2.40e-3},{6,9.07e-3},{7,2.38e-2},
    {8,5.00e-2},{9,8.96e-2},{10,1.44e-1},{11,2.14e-1},{12,3.00e-1},{13,4.00e-1},{14,5.14e-1},
    {15,6.41e-1},{16,7.81e-1},{17,9.31e-1},{18,1.09},{19,1.26},{20,1.44},{21,1.62},{22,1.82},
    {23,2.01},{24,2.22},{25,2.43},{26,2.64},{27,2.86},{28,3.08},{29,3.31},{30,3.54},{35,4.7},
    {40,6.0},{45,7.2},{50,8.5},{55,9.9}
  };
  const double _tol = std::pow(2, -53);
  double _mu, _t;
  sp_mat _A; 
  mat _B;
  int _s, _m_star;
  std::vector<std::vector<mat>> _Bs;

  //public:
  const int dim;
  mat result;

  SparseMatrixExponentialMultiply (
      const sp_mat& A, 
      const mat& B, 
      const double& t) 
    : dim(A.n_rows)
  {
    if (A.n_rows != A.n_cols) 
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Operator must be square");
    }

    if (A.n_cols != B.n_rows) 
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Dimension mismatch between operator and rhs");
    }

    _mu = trace(A) / double(A.n_rows);
    _A = A - _mu * speye(size(A)); 
    _B = B;
    _t = t;
    double _A_1_norm = norm(_A, 1);
    if (_t*_A_1_norm == 0.)
    {
      _m_star = 0;
      _s = 1;
    } else {
      OperatorNorm operator_norm(_t*_A, _t*_A_1_norm, 2);
      ivec frag = _fragment_3_1(operator_norm, _B.n_cols, _tol, 2);
      _m_star = frag[0];
      _s = frag[1];
      if (_s < 1 || _m_star < 1)
      {
        Rcpp::stop(prefix + "Matrix exponential failed");
      }
    }
    result = _expm_multiply_simple_core(_A, _B, _t, _mu, _m_star, _s, _tol, _Bs);
  }

  sp_mat reverse_differentiate (mat& dB, const mat& dX)
  {
    if (dX.n_rows != _B.n_rows || dX.n_cols != _B.n_cols)  
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Gradient dimensions do not match rhs");
    }
    return _expm_multiply_simple_core_reverse_differentiate(dX, dB, _A, _Bs, _mu, _s, _t);
  }

  sp_mat operator() (void)
  {
    return _A + _mu * speye(size(_A));
  }
  
  struct OperatorNorm
  {
    const sp_mat _A;
    double _A_1_norm, _scale;
    int _ell;
    std::map<int, double> _d;
    const bool _approximate_norm = true; //TODO: should allow this to be toggled
  
    OperatorNorm (const sp_mat& A, double A_1_norm, int ell, double scale = 1.)
      : _A (A), _A_1_norm (A_1_norm), _ell (ell), _scale (scale)
    {}

    double norm_1_power (const sp_mat& A, int p)
    {
      OneNormEst one_norm (A, p, _approximate_norm);
      return one_norm._est;
    }
  
    double norm_1 (void)
    {
      return _scale * _A_1_norm;
    }
  
    double d (int p)
    {
      if (_d.find(p) == _d.end())
      {
        double est = norm_1_power(_A, p);
        _d[p] = std::pow(est, 1.0 / double(p));
      }
      return _scale * _d[p];
    }
  
    double alpha (int p)
    {
      double a = d(p);
      double b = d(p+1);
      return a > b ? a : b;
    }
  };
  
  long double _compute_cost_div_m (int m, int p, OperatorNorm& operator_norm)
  {
    return (long double)(std::ceil(operator_norm.alpha(p) / _theta[m]));
  }
  
  int _compute_p_max (int m_max)
  {
    double sqrt_m_max = std::sqrt(m_max);
    int p_low = int(std::floor(sqrt_m_max));
    int p_high = int(std::ceil(sqrt_m_max + 1));
    int p_max = -1; 
    for (int p=p_low; p<p_high+1; ++p)
    {
      if (p*(p-1) <= m_max + 1)
      {
        if (p > p_max)
        {
          p_max = p;
        }
      }
    }
    return p_max;
  }
  
  bool _condition_3_13 (double A_1_norm, int n0, int m_max, int ell)
  {
    int p_max = _compute_p_max(m_max);
    int a = 2 * ell * p_max * (p_max + 3);
    double b = _theta[m_max] / double(n0 * m_max);
    return A_1_norm <= a * b;
  }
  
  ivec _fragment_3_1 (OperatorNorm& operator_norm, int n0, double tol, int ell)
  {
    if (ell < 1) 
    {
      Rcpp::stop(prefix + "ell < 1");
    }
    int m_max  = 55;
    long int best_m = -1;
    long double best_s = -1;
    if (_condition_3_13(operator_norm.norm_1(), n0, m_max, ell))
    {
      for (const auto& _theta_ptr : _theta)
      {
        long int m = _theta_ptr.first;
        double theta = _theta_ptr.second;
        long double s = std::ceil(operator_norm.norm_1() / theta);
        if (best_m == -1 || m*s < best_m*best_s)
        {
          best_m = m;
          best_s = s;
        }
      }
    } else {
      for (int p=2; p<_compute_p_max(m_max)+1; ++p)
      {
        for (int m=p*(p-1)-1; m<m_max+1; ++m)
        {
          if (_theta.find(m) != _theta.end())
          {
            long double s = _compute_cost_div_m(m, p, operator_norm);

            if (best_m == -1 || m*s < best_m*best_s)
            {
              best_m = m;
              best_s = s;
            }
          }
        }
      }
      best_s = std::max(best_s, (long double)(1.0));
    }
    return ivec({int(best_m), int(best_s)});
  }
  
  mat _expm_multiply_simple_core (
      const sp_mat& A, 
      mat B, 
      double t, 
      double mu, 
      int m_star, 
      int s, 
      double tol, 
      std::vector<std::vector<mat>>& Bs)
  {
    Bs.clear();
    Bs.resize(s+1);
    mat F = B; 
    double eta = std::exp(t * mu / double(s)); 
    for(int i=0; i<s; ++i)
    {
      Bs[i].reserve(m_star);
      B = F; 
      double c1 = norm(B, "inf");;
      for(int j=0; j<m_star; ++j)
      {
        double coeff = t / double(s*(j+1));
        Bs[i].push_back(B);
        B = coeff * A * B;
        double c2 = norm(B, "inf");
        F += B; 
        if (c1 + c2 <= tol * norm(F, "inf")) 
        {
          break;
        }
        c1 = c2;
      }
      F *= eta; 
    }
    Bs[s].push_back(F);
    return F; 
  }
  
  sp_mat _expm_multiply_simple_core_reverse_differentiate (
      const mat& dX, 
      mat& dB_out, 
      const sp_mat& A, 
      const std::vector<std::vector<mat>>& B, 
      const double mu, 
      const int s, 
      const double t)
  {
    double eta = std::exp(t * mu / double(s)); 
    double deta = 0;
    mat F = B[s][0];
    mat dF = dX;
    sp_mat dA = spones(A); 
    for(int i=s-1; i>=0; --i)
    {
      F /= eta;
      deta += dot(dF, F);
      dF = dF * eta;
      mat dB = zeros(size(dX));
      for(int j=B[i].size()-1; j>=0; --j)
      {
        // storing Bij seems costly but I don't think there's a way around it
        // as the only other option is repeatedly inverting A
        mat Bij = B[i][j]; //computed in forward pass; B at beginning of forward loop
        long double coeff = t / (long double)(s*(j+1));
        dB += dF;
        F = F - coeff * A * Bij; 
        for (sp_mat::iterator dA_ = dA.begin(); dA_ != dA.end(); ++dA_)
        {
          double Av = (*dA_) + coeff * 
            arma::dot(dB.row(dA_.row()), Bij.row(dA_.col()));
          if (Av == 0.0) 
          {
            Rcpp::stop(prefix + " Can't void nonzero element");
          }
          (*dA_) = Av;
        }
        dB = coeff * trans(A) * dB;
      }
      dF += dB;
    }
    //dmu = deta * eta * t/double(s);
    dB_out = dF; //note that it is NOT incremented!!!
    dA -= spones(A); 
    return dA;
  }
};

// ---------------------------------------------------

struct SparseMatrixExponentialMultiplySafe
{
  /*
   *  Split matrix exponential multiply into intervals, adaptively choosing
   *  the length of the interval. Assumes that the output should be left-stochastic.
   */

  // TODO:
  //  - Do we need to renormalise at end so that output columns sum to 1?

  const std::string prefix = "[SparseMatrixExponentialMultiplySafe] ";
  const double step_reduction = std::sqrt(0.5); // should this be settable?

  std::vector<SparseMatrixExponentialMultiply> mat_exp_mult;
  mat result;
  unsigned evals;

  SparseMatrixExponentialMultiplySafe (
      const sp_mat& A, 
      const mat& B, 
      const double& t,
      const double tol = 1e-12, //gives equivalent or better accuracy than expm::expm
      double step_min = 1.0e-4) 
  {
    if (t <= 0.0 || !std::isfinite(t))
    {
      Rcpp::stop(prefix + " Time step must be positive and finite");
    }

    if (step_min <= 0.0 || step_min > 1.0) 
    {
      Rcpp::stop(prefix + " Minimum step length must be in (0, 1]");
    }

    if (step_reduction <= 0.0 || step_reduction > 1.0) 
    {
      Rcpp::stop(prefix + " Step reduction factor must be in (0, 1]");
    }

    if (tol <= 0.0)
    {
      Rcpp::stop(prefix + " Convergence tolerance must be positive");
    }

    double step = t;
    double total = 0;
    unsigned n_steps = 0;
    mat _B = B;
    evals = 0;
    do 
    {
      double t_step = std::min(t - total, step);

      mat_exp_mult.emplace_back(A, _B, t_step);
      auto last_eval = mat_exp_mult.back();
      evals++;

      if (_check_left_stochastic(last_eval.result, tol))
      {
        _B = last_eval.result;
        total += t_step;
        if (total > t) Rcpp::stop(prefix + "Interval exceeded");
      } else {
        mat_exp_mult.pop_back();
        step *= step_reduction;
        if (step < step_min * t) 
        {
          A.save("_tempFile_coaldecoder_A_" + std::to_string(t) + ".txt", arma::coord_ascii); //DEBUG
          B.save("_tempFile_coaldecoder_B_" + std::to_string(t) + ".txt", arma::coord_ascii); //DEBUG
          Rcpp::stop(prefix + "Minimum allowed step length reached");
        }
      } 
    } while (total != t);

    auto last_eval = mat_exp_mult.back();
    result = last_eval.result;
  }

  bool _check_left_stochastic (const mat& x, const double& tol)
  {
    return all(vectorise(x) >= 0.0) && all(abs(sum(x, 0) - 1.0) <= tol);
  }

  sp_mat reverse_differentiate (mat& dB, const mat& dX)
  {
    sp_mat dA = spones(mat_exp_mult[0]._A);
    dB = dX;
    for (int i=mat_exp_mult.size()-1; i>=0; --i)
    {
      mat tmp_B;
      dA += mat_exp_mult[i].reverse_differentiate(tmp_B, dB);
      dB = tmp_B;
    }
    return dA - spones(mat_exp_mult[0]._A);
  }
};

//------------------- for debugging purposes

struct SparseMatrixExponentialMultiplyRescale
{
  //private:
  std::map<int, double> _theta = {
    {1,2.29e-16},{2,2.58e-8},{3,1.39e-5},{4,3.40e-4},{5,2.40e-3},{6,9.07e-3},{7,2.38e-2},
    {8,5.00e-2},{9,8.96e-2},{10,1.44e-1},{11,2.14e-1},{12,3.00e-1},{13,4.00e-1},{14,5.14e-1},
    {15,6.41e-1},{16,7.81e-1},{17,9.31e-1},{18,1.09},{19,1.26},{20,1.44},{21,1.62},{22,1.82},
    {23,2.01},{24,2.22},{25,2.43},{26,2.64},{27,2.86},{28,3.08},{29,3.31},{30,3.54},{35,4.7},
    {40,6.0},{45,7.2},{50,8.5},{55,9.9}
  };
  const double _tol = std::pow(2, -53);
  double _mu, _t;
  sp_mat _A; 
  mat _B;
  int _s, _m_star;
  std::vector<std::vector<mat>> _Bs;

  //public:
  const int dim;
  mat result;

  SparseMatrixExponentialMultiplyRescale (
      const sp_mat& A, 
      const mat& B, 
      const double& t)
    : dim(A.n_rows)
  {
    if (A.n_rows != A.n_cols) 
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Operator must be square");
    }

    if (A.n_cols != B.n_rows) 
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Dimension mismatch between operator and rhs");
    }

    _mu = trace(A) / double(A.n_rows);
    _A = A - _mu * speye(size(A)); 
    _B = B;
    _t = t;
    double _A_1_norm = norm(_A, 1);
    if (_t*_A_1_norm == 0.)
    {
      _m_star = 0;
      _s = 1;
    } else {
      OperatorNorm operator_norm(_t*_A, _t*_A_1_norm, 2);
      ivec frag = _fragment_3_1(operator_norm, _B.n_cols, _tol, 2);
      _m_star = frag[0];
      _s = frag[1];
    }
    // if _m_star / s doesn't fall within desired threshold, don't do anything, set flag
    // if (_m_star * s < )
    // valid = false;
    // um not much benefit to this nm
    result = _expm_multiply_simple_core(_A, _B, _t, _mu, _m_star, _s, _tol, _Bs);
  }

  sp_mat reverse_differentiate (mat& dB, const mat& dX)
  {
    if (dX.n_rows != _B.n_rows || dX.n_cols != _B.n_cols)  
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Gradient dimensions do not match rhs");
    }
    return _expm_multiply_simple_core_reverse_differentiate(dX, dB, _A, _Bs, _mu, _s, _t);
  }

  sp_mat operator() (void)
  {
    return _A + _mu * speye(size(_A));
  }

  double OpNormD (const int p) 
  {
    double _A_1_norm = norm(_A, 1);
    OperatorNorm operator_norm(_t*_A, _t*_A_1_norm, 2);
    return operator_norm.d(p);
  }
  
  struct OperatorNorm
  {
    const sp_mat _A;
    double _A_1_norm, _scale;
    int _ell;
    std::map<int, double> _d;
  
    OperatorNorm (const sp_mat& A, double A_1_norm, int ell, double scale = 1.)
      : _A (A), _A_1_norm (A_1_norm), _ell (ell), _scale (scale)
    {}

    double norm_1_power (const sp_mat& A, int p)
    {
      OneNormEst one_norm (A, p, true);
      return one_norm._est;
    }
  
    double norm_1 (void)
    {
      return _scale * _A_1_norm;
    }
  
    double d (int p)
    {
      if (_d.find(p) == _d.end())
      {
        double est = norm_1_power(_A, p);
        _d[p] = std::pow(est, 1.0 / double(p));
      }
      return _scale * _d[p];
    }
  
    double alpha (int p)
    {
      double a = d(p);
      double b = d(p+1);
      return a > b ? a : b;
    }
  };
  
  long double _compute_cost_div_m (int m, int p, OperatorNorm& operator_norm)
  {
    std::cout << "_compute_cost_div_m" << std::endl;//DEBUG
    std::cout << operator_norm.alpha(p) << " " << _theta[m] << std::endl;//DEBUG
    std::cout << std::ceil(operator_norm.alpha(p) / _theta[m]) << std::endl;//DEBUG
    std::cout << "/_compute_cost_div_m" << std::endl;//DEBUG
    return (long double)(std::ceil(operator_norm.alpha(p) / _theta[m]));
  }
  
  int _compute_p_max (int m_max)
  {
    double sqrt_m_max = std::sqrt(m_max);
    int p_low = int(std::floor(sqrt_m_max));
    int p_high = int(std::ceil(sqrt_m_max + 1));
    int p_max = -1; 
    for (int p=p_low; p<p_high+1; ++p)
    {
      if (p*(p-1) <= m_max + 1)
      {
        if (p > p_max)
        {
          p_max = p;
        }
      }
    }
    return p_max;
  }
  
  bool _condition_3_13 (double A_1_norm, int n0, int m_max, int ell)
  {
    int p_max = _compute_p_max(m_max);
    int a = 2 * ell * p_max * (p_max + 3);
    double b = _theta[m_max] / double(n0 * m_max);
    return A_1_norm <= a * b;
  }
  
  ivec _fragment_3_1 (OperatorNorm& operator_norm, int n0, double tol, int ell)
  {
    if (ell < 1) 
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] ell < 1");
    }
    int m_max  = 55;
    long int best_m = -1;
    long double best_s = -1;
    if (_condition_3_13(operator_norm.norm_1(), n0, m_max, ell))
    {
      for (const auto& _theta_ptr : _theta)
      {
        long int m = _theta_ptr.first;
        double theta = _theta_ptr.second;
        long double s = std::ceil(operator_norm.norm_1() / theta);
        if (best_m == -1 || m*s < best_m*best_s)
        {
          best_m = m;
          best_s = s;
        }
      }
    } else {
      for (int p=2; p<_compute_p_max(m_max)+1; ++p)
      {
        for (int m=p*(p-1)-1; m<m_max+1; ++m)
        {
          if (_theta.find(m) != _theta.end())
          {
            long double s = _compute_cost_div_m(m, p, operator_norm);

            if (best_m == -1 || m*s < best_m*best_s)
            {
              best_m = m;
              best_s = s;
            }
          }
        }
      }
      best_s = std::max(best_s, (long double)(1.0));
    }
    return ivec({int(best_m), int(best_s)});
  }
  
  mat _expm_multiply_simple_core (
      const sp_mat& A, 
      mat B, 
      double t, 
      double mu, 
      int m_star, 
      int s, 
      double tol, 
      std::vector<std::vector<mat>>& Bs)
  {
    Bs.clear();
    Bs.resize(s+1);
    mat F = B; 
    double eta = std::exp(t * mu / double(s)); 
    for(int i=0; i<s; ++i)
    {
      Bs[i].reserve(m_star);
      B = F; 
      double c1 = norm(B, "inf");;
      for(int j=0; j<m_star; ++j)
      {
        double coeff = t / double(s*(j+1));
        Bs[i].push_back(B);
        B = coeff * A * B;
        double c2 = norm(B, "inf");
        F += B; 
        if (c1 + c2 <= tol * norm(F, "inf")) 
        {
          break;
        }
        c1 = c2;
      }
      F *= eta; 
    }
    Bs[s].push_back(F);
    return F; 
  }
  
  sp_mat _expm_multiply_simple_core_reverse_differentiate (
      const mat& dX, 
      mat& dB_out, 
      const sp_mat& A, 
      const std::vector<std::vector<mat>>& B, 
      const double mu, 
      const int s, 
      const double t)
  {
    double eta = std::exp(t * mu / double(s)); 
    double deta = 0;
    mat F = B[s][0];
    mat dF = dX;
    sp_mat dA = spones(A); 
    for(int i=s-1; i>=0; --i)
    {
      F /= eta;
      deta += dot(dF, F);
      dF = dF * eta;
      mat dB = zeros(size(dX));
      for(int j=B[i].size()-1; j>=0; --j)
      {
        // I do think this could be optimized so that it only
        // needs to store the second index?
        // storing Bij seems costly but I don't think there's a way around it
        // as the only other option is repeatedly inverting A
        mat Bij = B[i][j]; //computed in forward pass; B at beginning of forward loop
        long double coeff = t / (long double)(s*(j+1));
        dB += dF;
        F = F - coeff * A * Bij; 
        for (sp_mat::iterator dA_ = dA.begin(); dA_ != dA.end(); ++dA_)
        {
          double Av = (*dA_) + coeff * 
            arma::dot(dB.row(dA_.row()), Bij.row(dA_.col()));
          if (Av == 0.0) 
          {
            Rcpp::stop("[SparseMatrixExponentialMultiply] Can't void nonzero element");
          }
          (*dA_) = Av;
        }
        dB = coeff * trans(A) * dB;
      }
      dF += dB;
    }
    //dmu = deta * eta * t/double(s);
    dB_out = dF; //note that it is NOT incremented!!!
    dA -= spones(A); 
    return dA;
  }
};

#endif
