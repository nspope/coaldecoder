#ifndef _EXPM_REVDIF_SPARSE
#define _EXPM_REVDIF_SPARSE

#include <RcppArmadillo.h> 
#include <map>

// C++ implementation based on:
// https://github.com/scipy/scipy/blob/v1.8.0/scipy/sparse/linalg/_expm_multiply.py#L56-L142

// reverse diff is by Nate Pope

using namespace arma;

struct SparseMatrixExponentialMultiply
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
    }
    result = _expm_multiply_simple_core(_A, _B, _t, _mu, _m_star, _s, _tol, _Bs);
  }

  sp_mat reverse_differentiate (mat& dB, const mat& dX)
  {
    if (dX.n_rows != _B.n_rows || dX.n_cols != _B.n_cols)  
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Gradient dimensions do not match rhs");
    }
    if (dB.n_rows != _B.n_rows || dB.n_cols != _B.n_cols)  
    {
      Rcpp::stop("[SparseMatrixExponentialMultiply] Output dimensions do not match rhs");
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
  
    OperatorNorm (const sp_mat& A, double A_1_norm, int ell, double scale = 1.)
      : _A (A), _A_1_norm (A_1_norm), _ell (ell), _scale (scale)
    {}

    double norm_1_power (const sp_mat& A, int p)
    {
      rowvec x = ones<rowvec>(A.n_rows);
      for (int i=0; i<p; ++i)
      {
        x *= A;
      }
      return abs(x).max();
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
  
  long int _compute_cost_div_m (int m, int p, OperatorNorm& operator_norm)
  {
    return (long int)(std::ceil(operator_norm.alpha(p) / _theta[m]));
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
    long int best_s = -1;
    if (_condition_3_13(operator_norm.norm_1(), n0, m_max, ell))
    {
      for (const auto& _theta_ptr : _theta)
      {
        long int m = _theta_ptr.first;
        double theta = _theta_ptr.second;
        long int s = std::ceil(operator_norm.norm_1() / theta);
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
            long int s = _compute_cost_div_m(m, p, operator_norm);

            if (best_m == -1 || m*s < best_m*best_s)
            {
              best_m = m;
              best_s = s;
            }
          }
        }
      }
      best_s = std::max(best_s, (long int)(1));
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
    dB_out += dF; //note that it is incremented!!!
    dA -= spones(A); 
    return dA;
  }
};

#endif
