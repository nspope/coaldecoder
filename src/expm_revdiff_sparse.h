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

// C++ implementation based on:
// https://github.com/scipy/scipy/blob/v1.8.0/scipy/sparse/linalg/_expm_multiply.py#L56-L142

// reverse diff is by Nate Pope

using namespace arma;

struct OneNormEst
{
  /*
   *  Calculate 1-norm of power of sparse matrix, exactly or approximately
   */

  const std::string prefix = "[OneNormEst] ";
  const double eps = 4 * datum::eps;
  const unsigned iter_max = 10;

  ivec _v;
  vec _w;
  uword _iter;
  double _est;

  OneNormEst (const sp_mat& A, const int p=0, const bool verbose=false, const bool approximate=true)  
  {
    //compute exact if !approximate or condition met?
    _est = approximate ? approximate_norm(A, p, verbose) : exact_norm(A, p, verbose);
  }

  double exact_norm (const sp_mat& A, const int p, const bool verbose)
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

  double approximate_norm (const sp_mat& A, const int p, const bool verbose)
  {
    uword n = A.n_rows;
    uword t = std::min(int(n), 5);

    if (A.n_rows != A.n_cols)
    {
      Rcpp::stop(prefix + "Matrix must be square");
    }

    if (t < 1 || iter_max < 1) Rcpp::stop(prefix + "Invalid inputs");

    mat X = randu(n, t); 
    X.each_col([](vec& x) { x /= accu(x); });
    uvec been_there = zeros<uvec>(n);
    mat I_t = eye(t, t);
    double est_old = 0.0;
    double est = 0.0;
    mat S = zeros(n, t);
    mat S_old = S;
    vec w = zeros(n);
    uword imax;
    unsigned iter;
    for (iter=0; iter<iter_max+1; ++iter)
    {
      mat Y = A_x(A, X, p);
      rowvec cY = sum(abs(Y), 0);
      imax = cY.index_max();
      est = cY[imax];
      if (est > est_old || iter == 1)
      {
        w = Y.col(imax);
      }
      if (iter >= 1 && est < est_old)
      {
        est = est_old;
        break;
      }
      est_old = est;
      S_old = S;
      if (iter == iter_max)
      {
        if (verbose)
        {
          Rcpp::Rcout << prefix << "Did not converge in " << iter_max << " iterations" << std::endl;
        }
        break;
      }
      S = sign(Y);
      mat S_check = abs(trans(S_old) * S - n);
      uvec partest (S_check.n_cols);
      for (uword i=0; i<S_check.n_cols; ++i)
      {
        uvec check = find(S_check.col(i) < eps * n);
        partest.at(i) = check.n_elem;
      }
      if (all(partest > 0))
      {
        if (verbose)
        {
          Rcpp::Rcout << prefix << "Hit a cycle (1), stopping iterations" << std::endl;
        }
        break;
      }
      if (any(partest > 0))
      {
        uvec check = find(partest > 0);
        for (auto j : check)
        {
          S.col(j) = sign(randu(n) - 0.5);
        }
      }
      S_check = trans(S) * S - I_t;
      for (uword i=0; i<S_check.n_cols; ++i)
      {
        uvec check = find(S_check.col(i) == n);
        partest.at(i) = check.n_elem;
      }
      if (any(partest > 0))
      {
        uvec check = find(partest > 0);
        for (auto j : check)
        {
          S.col(j) = sign(randu(n) - 0.5);
        }
      }
      mat Z = At_x(A, S, p);
      mat h = abs(Z);
      h.transform([](double x){ return(std::max(2.0, x)); });
      uvec mhi (h.n_cols);
      for (unsigned i=0; i<h.n_cols; ++i)
      {
        mhi[i] = h.col(i).index_max();
      }
      if (iter >= 1 && all(mhi == imax))
      {
        if (verbose)
        {
          Rcpp::Rcout << prefix << "Hit a cycle (2), stopping iterations" << std::endl;
        }
        break;
      }
      umat indmat (size(h));
      for (unsigned i=0; i<h.n_cols; ++i)
      {
        indmat.col(i) = sort_index(h.col(i), "descend");
        h.col(i) = sort(h.col(i), "descend");
      }
      uvec ind = vectorise(indmat);
      if (t > 1)
      {
        uvec firstind = ind.head(t);
        if (all(been_there.elem(firstind)))
        {
          break;
        }
        arma::uvec keep_ind = find(been_there.elem(ind) == 0);
        ind = ind.elem(keep_ind);
        if (ind.n_elem < t)
        {
          if (verbose)
          {
            Rcpp::Rcout << prefix << "Not enough new vectors, stopping iterations" << std::endl;
          }
          break;
        }
      }
      X = zeros(n, t);
      for (unsigned i=0; i<t; ++i)
      {
        X.at(ind[i], i) = 1.0;
        been_there[ind[i]] = 1;
      }
    }

    ivec v = zeros<ivec>(n);
    v[imax] = 1;

    _v = v;
    _w = w;
    _iter = iter;

    return est;
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
    const bool _verbose = false;
  
    OperatorNorm (const sp_mat& A, double A_1_norm, int ell, double scale = 1.)
      : _A (A), _A_1_norm (A_1_norm), _ell (ell), _scale (scale)
    {}

    double norm_1_power (const sp_mat& A, int p)
    {
      OneNormEst one_norm (A, p, _verbose, true);
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
          Rcpp::stop(prefix + "Minimum allowed step length reached");
        }
      } 
    } while (total != t);

    auto last_eval = mat_exp_mult.back();
    result = last_eval.result;
  }

  bool _check_left_stochastic (const mat& x, const double& tol)
  {
    return all(vectorise(x) >= 0) && all(abs(sum(x, 0) - 1.0) <= tol);
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

//-------------------
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

  double OpNormD (const int p) //DEBUG
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
      //rowvec x = ones<rowvec>(A.n_rows);
      //for (int i=0; i<p; ++i)
      //{
      //  x *= A;
      //}
      //return abs(x).max();
      OneNormEst one_norm (A, p, true, true);
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
