#ifndef _EXPM_REVDIF
#define _EXPM_REVDIF

#include <RcppArmadillo.h> 
#include <map>

// C++ implementation based on:
// https://github.com/scipy/scipy/blob/v1.8.0/scipy/sparse/linalg/_expm_multiply.py#L56-L142

// reverse diff is by Nate Pope

using namespace arma;

struct matrix_exponential_multiply
{
  public:
  const int dim;
  mat result;

  matrix_exponential_multiply (const mat& A, const mat& B, const double& t) : dim(A.n_rows)
  {
    if (A.n_rows != A.n_cols) Rcpp::stop("A must be square");
    if (A.n_cols != B.n_rows) Rcpp::stop("dimension mismatch between A and B");
    _mu = _trace(A) / double(A.n_rows);
    _A = A - _mu * eye(size(A));
    _B = B;
    _t = t;
    double _A_1_norm = _exact_1_norm(_A);
    if (_t*_A_1_norm == 0.)
    {
      _m_star = 0;
      _s = 1;
    } else {
      _operator_norm_info norm_info(_t*_A, _t*_A_1_norm, 2);
      ivec frag = _fragment_3_1(norm_info, _B.n_cols, _tol, 2);
      _m_star = frag[0];
      _s = frag[1];
    }
    result = _expm_multiply_simple_core(_A, _B, _t, _mu, _m_star, _s, _tol, _B_end);
  }

  void reverse_differentiate (mat& gradient_wrt_A, mat& gradient_wrt_B, const mat& gradient)
  {
    if (gradient.n_rows != _B.n_rows || gradient.n_cols != _B.n_cols) Rcpp::stop("gradient does not match input");
    double dmu = 0; //evidently cancels, but keeping just in case
    _expm_multiply_simple_core_revdif(gradient, gradient_wrt_A, gradient_wrt_B, dmu, _A, _B_end, _mu, _s, _t);
  }

  mat operator() (void)
  {
    return _A + _mu * eye(size(_A));
  }

  private:
  std::map<int, double> _theta = {
    {1,2.29e-16},{2,2.58e-8},{3,1.39e-5},{4,3.40e-4},{5,2.40e-3},{6,9.07e-3},{7,2.38e-2},
    {8,5.00e-2},{9,8.96e-2},{10,1.44e-1},{11,2.14e-1},{12,3.00e-1},{13,4.00e-1},{14,5.14e-1},
    {15,6.41e-1},{16,7.81e-1},{17,9.31e-1},{18,1.09},{19,1.26},{20,1.44},{21,1.62},{22,1.82},
    {23,2.01},{24,2.22},{25,2.43},{26,2.64},{27,2.86},{28,3.08},{29,3.31},{30,3.54},{35,4.7},
    {40,6.0},{45,7.2},{50,8.5},{55,9.9}
  };
  const double _tol = std::pow(2, -53);
  double _mu, _t;
  mat _A, _B;
  int _s, _m_star;
  std::vector<std::vector<mat>> _B_end;
  
  double _exact_1_norm (const mat& A)
  {
    return norm(A, 1);
  }
  
  double _trace (const mat& A)
  {
    return trace(A);
  }
  
  double _exact_inf_norm (const mat& A)
  {
    return norm(A, "inf");
  }
  
  struct _operator_norm_info
  {
    mat _A;
    double _A_1_norm, _scale;
    int _ell;
    std::map<int, double> _d;
  
    _operator_norm_info (mat A, double A_1_norm, int ell, double scale = 1.)
      : _A (A), _A_1_norm (A_1_norm), _ell (ell), _scale (scale)
    {}

// here is rough conversion of onenormest, need this when A gets huge
//def onenormest(A, t=2, itmax=5, compute_v=False, compute_w=False):
//
//  if (A.n_rows != A.n_cols) Rcpp::stop("A must be square");
//  int n = A.n_cols;
//  vec v, w;
//  double est;
//  if (t >= n) {
//    rowvec col_abs_sums = sum(abs(A), 0);
//    int argmax_j = col_abs_sums.index_max();
//    v = zeros(n);
//    v[argmax_j] = 1.;
//    w = A.col(j);
//    est = col_abs_sums[j];
//  } else {
//    _onenormest_core(est, v, w, nmults, nresamples, A, trans(A), t, itmax);
//  }
//  return est;
//
//def _onenormest_core(A, AT, t, itmax):
//  if (itmax < 2) Rcpp::stop("at least two iterations are required");
//  if (t < 1) Rcpp::stop("at least one column is required");
//  int n = A.n_rows;
//  if (t >= n) Rcpp::stop("t should be smaller than the order of A");
//  int nmults = 0, nresamples = 0;
//  vec X = ones(n, t);
//  if (t > 1)
//  {
//    for (int i=1; i<t; ++i)
//    {
//      resample_column(i, X)
//    }
//    for (int i=0; i<t; ++i)
//    {
//      while (column_needs_resampling(i, X))
//      {
//        resample_column(i, X);
//        nresamples++;
//      }
//    }
//    X /= double(n);
//    ivec ind_hist = zeros(0); //ind_hist = np.zeros(0, dtype=np.intp)
//    double est_old = 0.;
//    mat S = zeros(n, t);
//    int k = 1;
//    int ind = -1; //ind = None ... hmmm
//    int best_j;
//    vec w, h;
//    while (true) {
//      mat Y = A * X;
//      nmults++;
//      double mags = _sum_abs_axis0(Y);
//      est = mags.max();
//      best_j = mags.index_max();
//      if (est > est_old || k == 2)
//      {
//        if (k >= 2)
//        {
//          ind_best = ind.at(best_j);
//        }
//        w = Y.col(best_j);
//      }
//      if (k >= 2 && est <= est_old)
//      {
//        est = est_old;
//        break;
//      }
//      est_old = est;
//      mat S_old = S;
//      if (k > itmax) break;
//      S = sign_round_up(Y); //also had del Y
//      if (every_col_of_X_is_parallel_to_a_col_of_Y(S, S_old))
//        break;
//      if (t > 1)
//      {
//        for (int i=0; i<t; ++i)
//        {
//          while (column_needs_resampling(i, S, S_old))
//          {
//            resample_column(i, S);
//            nresamples++;
//          }
//        }
//      }
//      //del S_old
//      mat Z = arma::trans(A) * S;
//      nmults++;
//      h = _max_abs_axis1(Z);
//      //del Z
//      if (k >= 2 && h.max() == h.at(ind_best)) break;
//      ind = sort_index(h, "descend");
//      ind = ind.head(t+ind_hist.n_elem);
//      //ind = np.argsort(h)[::-1][:t+len(ind_hist)].copy()
//      //del h
//      uvec ind_head_t = ind.head(t);
//      uvec seen_head_t = ind_head_t in ind_hist;
//      if (t > 1)
//      {
//        if (all(seen_head_t == 1)) break;
//        uvec seen = ind in ind_hist;
//        ind = join_vert(ind.elem(find(seen == 0)), 
//                        ind.elem(find(seen == 1)));
//      }
//      for (int j=0; j<t; ++j)
//      {
//        X.col(j).zeros();
//        X.at(ind.at(j),j) = 1.;
//      }
//      new_ind = ind_head_t.elem(find(seen_head_t == 0));
//      ind_hist = join_vert(ind_hist, new_ind);
//      k += 1;
//    }
//    v = zeros(n);
//    v.at(ind_best) = 1.;
//    return est, v, w, nmults, nresamples;

    double _onenormest_matrix_power(const mat& A, int p, int t=2)
    {
      //t is evidently unused
      return norm(powmat(A, p), 1);
    }
  
    void set_scale (double scale)
    {
      _scale = scale;
    }
  
    double onenorm (void)
    {
      return _scale * _A_1_norm;
    }
  
    double d (int p)
    {
      if (_d.find(p) == _d.end())
      {
        double est = _onenormest_matrix_power(_A, p, _ell);
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
  
  long int _compute_cost_div_m (int m, int p, _operator_norm_info& norm_info)
  {
    return (long int)(std::ceil(norm_info.alpha(p) / _theta[m]));
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
  
  ivec _fragment_3_1 (_operator_norm_info& norm_info, int n0, double tol, int ell)
  {
    if (ell < 1) Rcpp::stop("ERROR ell < 1");
    int m_max  = 55;
    long int best_m = -1;
    long int best_s = -1;
    if (_condition_3_13(norm_info.onenorm(), n0, m_max, ell))
    {
      for (const auto& _theta_ptr : _theta)
      {
        long int m = _theta_ptr.first;
        double theta = _theta_ptr.second;
        long int s = std::ceil(norm_info.onenorm() / theta);
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
            long int s = _compute_cost_div_m(m, p, norm_info);

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
  
  mat _expm_multiply_simple_core (const mat& A, mat B, double t, double mu, int m_star, int s, double tol, 
                                  std::vector<std::vector<mat>>& B_end)
  {
    B_end.resize(s+1);
    mat F = B; 
    double eta = std::exp(t * mu / double(s)); 
    for(int i=0; i<s; ++i)
    {
      B_end[i].reserve(m_star);
      B = F; 
      double c1 = _exact_inf_norm(B);
      for(int j=0; j<m_star; ++j)
      {
        double coeff = t / double(s*(j+1));
        B_end[i].push_back(B);
        B = coeff * A * B;
        double c2 = _exact_inf_norm(B);
        F += B; 
        if (c1 + c2 <= tol * _exact_inf_norm(F)) 
        {
          break;
        }
        c1 = c2;
      }
      F *= eta; 
    }
    B_end[s].push_back(F);
    return F; 
  }
  
  //mat _expm_multiply_simple (mat A, mat B, const double t, int& s, std::vector<std::vector<mat>>& B_end)
  //{
  //  mat ident = eye(size(A));
  //  int n = A.n_rows;
  //  int n0 = B.n_cols;
  //  double tol = std::pow(2, -53);
  //  double mu = _trace(A) / double(n);
  //  A -= mu * ident; 
  //  double A_1_norm = _exact_1_norm(A);
  //  int m_star;
  //  if (t*A_1_norm == 0.)
  //  {
  //    m_star = 0;
  //    s = 1;
  //  } else {
  //    int ell = 2;
  //    _operator_norm_info norm_info(t*A, t*A_1_norm, ell);
  //    ivec frag = _fragment_3_1(norm_info, n0, tol, ell);
  //    m_star = frag[0];
  //    s = frag[1];
  //  }
  //  return _expm_multiply_simple_core(A, B, t, mu, m_star, s, tol, B_end);
  //}
  
  void _expm_multiply_simple_core_revdif (const mat& dX, mat& dA, mat& dB, double& dmu, const mat& A, const std::vector<std::vector<mat>>& B, const double mu, const int s, const double t)
  {
    double eta = std::exp(t * mu / double(s)); 
    double deta = 0;
    mat F = B[s][0];
    mat dF = dX;
    dA = zeros(size(A));
    for(int i=s-1; i>=0; --i)
    {
      F /= eta;
      deta += dot(dF, F);
      dF = dF * eta;
      dB = zeros(size(dX));
      for(int j=B[i].size()-1; j>=0; --j)
      {
        mat Bij = B[i][j]; //computed in forward pass; B at beginning of forward loop
        long double coeff = t / (long double)(s*(j+1));
        dB += dF;
        F = F - coeff * A * Bij;
        dA += dB * trans(Bij) * coeff;
        dB = coeff * trans(A) * dB;
      }
      dF += dB;
    }
    dmu = deta * eta * t/double(s);
    dB = dF;
  
    //forward:
    //F[0,0] = in;
    //for(int i=0; i<s; ++i)
    //{
    //  B[i,0] = F[i,0];
    //  for(int j=0; j<=m_break[i]; ++j)
    //  {
    //    B[i,j+1] = coeff[i,j] * A * B[i,j]; 
    //    F[i,j+1] = F[i,j] + B[i,j+1]; 
    //  }
    //  F[i+1,0] = F[i,m_break[i]+1] * eta; 
    //}
    //out:F[s,0]
  
    //backward:
    // dF[s,0] = in;
    // for(int i=s-1; i>=0; --i)
    // {
    //   deta[i]  = dot(dF[i+1,0], F[i,m_break[i]+1]);
    //   dF[i,m_break[i]+1] = dF[i+1,0] * eta;
    //   dB[i,m_break[i]+1] = zeros;
    //   for(int j=m_break[i]; j>=0; --j)
    //   {
    //     dF[i,j] = dF[i,j+1];
    //     dB[i,j+1] += dF[i,j+1];
    //     dA += dB[i,j+1] * trans(B[i,j]) * coeff[i,j];
    //     dB[i,j] = coeff[i,j] * A.trans() * dB[i,j+1];
    //   }
    //   dF[i,0] += dB[i,0];
    // }
    // out:dF[0,0]
  }
  
  //arma::mat _expm_multiply_simple_revdif(const mat& dX, mat& dA, mat& dB, mat A, const std::vector<std::vector<mat>>& B, const double t, const int s) 
  //{
  //  int n = A.n_rows;
  //  double dmu = 0;
  //  double mu = _trace(A)/double(n); 
  //  mat ident = eye(size(A));
  //  A -= mu * ident;
  //  _expm_multiply_simple_core_revdif(dX, dA, dB, dmu, A, B, mu, s, t);
  //  //dmu -= accu(dA.diag());
  //  //dA.diag() -= dmu/double(n);
  //}
};

//Rcpp::List my_expm (arma::mat A, arma::mat B, double t)
//{
//  std::vector<std::vector<mat>> B_end;
//  int s;
//  mat out = _expm_multiply_simple (A, B, t, s, B_end);
//  mat dA; mat dB; mat dX = B;
//  _expm_multiply_simple_revdif (dX, dA, dB, A, B_end, t, s);
//  return Rcpp::List::create(
//      Rcpp::_["dA"] = dA,
//      Rcpp::_["dB"] = dB,
//      Rcpp::_["X"] = out,
//      Rcpp::_["s"] = s
//      );
//}

#endif
