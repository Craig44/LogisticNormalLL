// A TMB template script to calculate the negative loglikelihood
// and simulate from a Logistic Normal likelihood

// Forward Declaration of functions that are used in NLLlogistnorm and described below the NLLlogistnorm function
template <class Type> 
vector<Type>  RecursiveFilter(vector<Type> ar_coef, int nBins, vector<Type> initial_vals);
template <class Type> 
vector<Type> ARMAacf(Type AR, Type MA, int nBin);
template <class Type> 
vector<Type> ARMAacf_non_arma(vector<Type> AR, int nBin);
template <class Type> 
matrix<Type>  covariance_logistic(Type sigma, vector<Type> phi, int N_bins, int ARMA);
/*
 * The main functions which returns the Negative loglikelihood for the Logistic Normal according to Francis 2014
 */
template <class Type>  
Type NLLlogistnorm(array<Type>& obs_mat, matrix<Type>& exp_mat, vector<Type>& N,  Type sigma, vector<Type> phi, vector<Type> bin_labels, int ARMA) {
  int A = bin_labels.size();
  Type Y = Type(exp_mat.rows());
  int N_bins = exp_mat.cols();
  vector<Type> weights_by_year(exp_mat.rows());
  Type mean_N = N.mean();
  bool weights_the_same = true;
  // Calculate a weights that are essentially scaled N
  for (int this_N = 0; this_N < N.size();++this_N) {
    weights_by_year[this_N] = sqrt(mean_N/N[this_N]);
    if (this_N > 1) {
      if (weights_by_year[this_N] != weights_by_year[this_N - 1])
        weights_the_same = false;
    }
  }
  
  matrix<Type> logistic_covariance(N_bins,N_bins);
  logistic_covariance.setZero();
  matrix<Type> k_mat(N_bins - 1,N_bins);
  k_mat.setZero();
  matrix<Type> V_mat(N_bins - 1,N_bins);
  V_mat.setZero();
  matrix<Type> ww(N_bins,N_bins - 1);
  ww.setZero();
  //for (int a =0; a < N_bins; ++a)
  // logistic_covariance(a,a) = sigma * sigma;
  
  logistic_covariance = covariance_logistic(sigma,phi,N_bins,ARMA);
  // Get the covariance matrix with given parameters Sigma, phi and ARMA
  //logistic_covariance = covariance_logistic(sigma,phi,N_bins,ARMA);
  
  for (int i = 0; i < (N_bins - 1); ++i) {
    k_mat(i,i) = 1.0;
    k_mat(i, A - 1) = -1.0;
  }
  matrix<Type> t_kmat = logistic_covariance * k_mat.transpose();

  V_mat = k_mat * t_kmat;

  // Calculate log determinant and inverse
  Type log_det = atomic::logdet(V_mat);
  Type temp_val = 0;
  matrix<Type> V_mat_inv = atomic::matinvpd(V_mat,temp_val);
  Type log_obs_tot = obs_mat.log().sum();
  Type neg_ll_LN = 0.5 * Y * Type(A - 1) * log(2*3.1415926535) + log_obs_tot + 0.5 * Y * log_det;
  
  if (not weights_the_same) {
    neg_ll_LN += (Type(N_bins) - 1) * weights_by_year.log().sum();
  }
  
  
  for (int i = 0; i < exp_mat.rows(); ++i) {
    for (int j = 0; j < N_bins - 1; ++j) {
      Type l_obs = obs_mat(i,j) / obs_mat(i,N_bins - 1);
      Type l_exp = exp_mat(i,j) / exp_mat(i,N_bins - 1);
      ww(i, j) = log(l_obs) - log(l_exp);
    }
  }

  // TODO add robust which is what is in the C++ code
  Type temp2;
  matrix<Type> year_val, temp1;
  for (unsigned year = 0; year < exp_mat.rows(); ++year) {
    year_val = ww.row(year);
    temp1 = year_val * V_mat_inv * year_val.transpose();// * year_val;
    temp2 = 0.5 / (weights_by_year[year] * weights_by_year[year]);
    neg_ll_LN += (temp1.sum() * temp2);
  }
  return neg_ll_LN;
}

/**
 * @covariance_logistic - calculate the covariance matrix that is used in the logistic normal
 */
template <class Type> 
matrix<Type>  covariance_logistic(Type sigma, vector<Type> phi, int N_bins,int ARMA) {
  unsigned n_phi = phi.size();
  matrix<Type> covar(N_bins,N_bins);
  
  covar.setZero();
  if (phi[0] == 0.0 & (n_phi == 1)) {
    // no auto correlation
    for (unsigned diag = 0; diag < N_bins; ++ diag)
      covar(diag, diag) = sigma * sigma;
  } else {
    //covar = calculate_multivariate_normal_covariance(sigma, phi, N_bins, ARMA, false, bin_labels);
    vector<Type> rho(N_bins - 1);
    
    if (n_phi == 1) {
      for(int i = 1; i <= N_bins - 1; i++) {
        rho[i - 1]= pow(phi[0],(Type)i);
      }
    } else if (n_phi == 2 && (ARMA == 1)) {
      rho = ARMAacf(phi[0], phi[1], N_bins);
    } else {
      rho = ARMAacf_non_arma(phi, N_bins);
    }
    // Simple unsexed example currently
    // Create an identy matrix.
    for (int diag = 0; diag < N_bins; ++ diag) {
      covar(diag,diag) = 1.0 * sigma * sigma;
    }
    for (int diag = 0; diag < N_bins; ++ diag) {
      for (int row = 0; row < N_bins; ++ row) {
        for (int col = 0; col < N_bins; ++ col) {
          if (fabs(row - col) == diag + 1) {
            covar(row, col) = rho[diag] * sigma * sigma;
            covar(col, row) = rho[diag] * sigma * sigma; // make it symetric
            
          }
        }
      }
    }
  }
  return covar;
}

/**
 * @RecursiveFilter 
 * @param ar_coef an AR(2) coefficient
 * @param nBin number of ages
 * @param initial_vals initial vals
 */
template <class Type> 
vector<Type>  RecursiveFilter(vector<Type> ar_coef, int nBins, vector<Type> initial_vals) {
  vector<Type> store_vec(nBins + 1);
  
  if (ar_coef.size() > 2) {
    std::cerr <<  "RecursiveFilter(): has not been coded for more than 2 AR coeffiecients" << std::endl;
  }
  
  store_vec[0] = initial_vals[1];
  store_vec[1] = initial_vals[0];
  vector<Type> result(store_vec.size() - 1);
  for (unsigned i = 1; i < nBins + 1; ++i) {
    if (i == 1) {
      store_vec[i] =   store_vec[i - 1] *ar_coef[0]  + store_vec[i] *  ar_coef[1];
    } else {
      store_vec[i] = store_vec[i - 1] *  ar_coef[0] + store_vec[i - 2] * ar_coef[1];
    }
  }
  for (unsigned i = 1; i < store_vec.size(); ++i)
    result[i - 1] = store_vec[i];
  
  return result;
}

// Compute the theoretical autocorrelation function or partial autocorrelation function for an ARMA process.
// The methods used follow Brockwell & Davis (1991, section 3.3). Their equations (3.3.8) are solved for the
// autocovariances at lags 0, …, max(p, q+1), and the remaining autocorrelations are given by a recursive filter.
/**
* @description Calculate the auto correlation for the covariance matrix
* @param ar_coef an AR(2) coefficient
* @param ma_coef an MA(2) coefficien
* @param nBin number of ages
*/
template <class Type> 
vector<Type> ARMAacf(Type AR, Type MA, int nBin) {
  std::cout << "Beginning method ARMAacf() " << std::endl;
  int p = 1;
  int q = 1;
  vector<Type> MA_coef(1);
  MA_coef(0) = MA;
  
  if (!p && !q) {
    std::cerr << "empty model supplied" << std::endl;
  }
  int r = fmax(p, q + 1);
  
  vector<Type> AR_coef(1 + (r - p));
  AR_coef(0) = AR;
  
  for (int i = 0; i <= (r - p); ++i) {
    AR_coef(i + 1) = 0.0;
    p = r;
  }
  
  matrix<Type> A1(p + 1,2 * p + 1), ind(2 * p + 1,p + 1);
  A1.setZero();
  ind.setZero();
  for (int i = 0; i < ind.rows(); ++i) {
    for (int j = 0; j < ind.cols(); ++j) {
      ind(i,0) = Type(i + 1);
      ind(i, 1) = Type((i + 1) + (j + 1) - (i + 1));
    }
  }
  for (unsigned i = 0; i < A1.rows(); ++i) {
    A1(i ,i) = 1.0;
    A1(i, i + 1) = -AR_coef[0];
  }
  
  A1(0 ,A1.rows() - 1) = -AR_coef[1];
  A1(A1.rows() - 1, 0) = -AR_coef[1];
  A1(A1.rows() - 1, 1) = -AR_coef[0];
  
  vector<Type> psi(2);
  psi(0) = Type(1.0);
  psi(1) = AR + MA;
  
  vector<Type> theta(q + 3);
  theta(0) = Type(1.0);
  theta(1) = MA;
  for (int i = 0; i <= q; ++i)
    theta(i + 2) = Type(0.0);
  // Declare Eigen variables
  
  matrix<Type> A_eig(3,3); // 3f means 3x3 elements of type Type
  A_eig.setZero();
  matrix<Type> rhs(3,1);
  rhs.setZero();
  vector<int> seq(p + 1);
  
  // Calculate rhs
  
  for (int i = 0; i <= q; ++i) {
    Type x1 ,x2;
    x1 = psi[0]*theta[i];
    x2 = psi[1]*theta[i + q];
    rhs(i,0) = x1 + x2;
  }
  rhs(2,0) = 0.0;
  
  // Use the eigen library yo solve the inverse of for A with known vector B
  //vector<Type> Ind;
  for (int i = 0; i <= p; ++i) {
    seq(i) = (p - i);
  }
  
  for (int i = 0; i <= p; ++i) {
    for (int j = 0; j <= p; ++j) {
      A_eig(i,j) = Type(A1(seq[i], seq[j]));
    }
  }
  // Solve for A_eig given rhs
  //Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> x = A_eig.matrix().colPivHouseholderQr().solve(rhs);
  matrix<Type> iA = atomic::matinv(A_eig);
  // Populate A_eig and rhs
  matrix<Type> x = atomic::matmul(iA, rhs);

  vector<Type> final_acf(3);//,Cor;
  for (unsigned i = 1; i <= 2; ++i) {
    final_acf(i -1) = (x(i,0) / x(0,0));
  }
  
  //cout << "solution = " << x << endl;
  vector<Type> Cor = RecursiveFilter(AR_coef, nBin, final_acf);
  vector<Type> newVec(nBin);
  for (unsigned i = 0; i < nBin; ++i) {
    if (i == 0)
      newVec[i] = final_acf[0];
    else if (i == 1) {
      newVec[i] = final_acf[1];
    } else {
      newVec[i] = Cor[i];
    }
  }
  return newVec;
}

/**
 * @description An overloaded version of the above function, that only takes a single AR parameters.
 * @param AR an AR( vector of AR(1) coeffecients
 * @param nBin number of ages
 */
template <class Type> 
vector<Type> ARMAacf_non_arma(vector<Type> AR, int nBin) {
  unsigned p = AR.size();
  
  if (p > 2) {
    std::cerr << "This function has not been coded for more the two AR coeffs.";
  }
  int r = p;
  matrix<Type> A1(p + 1,2 * p + 1);
  A1.setZero();
  for (unsigned i = 0; i < A1.rows(); ++i) {
    A1(i ,i) = 1.0;
    A1(i, i + 1) = -AR[0];
    A1(i, i + 2) = -AR[1];
    
  }
  vector<Type> rhs(p + 1);
  rhs[0] = 1.0;
  
  vector<int> seq(p + 1);
  for (int i = 0; i <= p; ++i) {
    seq[i] = (p - i);
  }
  
  matrix<Type> A_temp(p + 1,p + 1);
  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p ; ++j) {
      if (j == 2)
        A_temp(i,j) = A1(i,j);
      else
        A_temp(i,j) = A1(i,j) + A1(i,2 * p  - j);
    }
  }
  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p ; ++j) {
      A_temp(i,j) =  A_temp(seq[i],seq[j]);
    }
  }
  // the bodge
  A_temp(1,2) = 0.0;
  Type log_det = 0.0;
  matrix<Type>A_inv = atomic::matinvpd(A_temp,log_det);
  
  
  // Take the first column of the inverted matrix
  vector<Type> final_acf(2);
  vector<Type> Acf(p + 1);
  
  for (unsigned i = 0; i < p + 1; ++i) {
    Acf(i) = A_inv(i,0);
  }
  
  // Divide the last two elements by the first element.
  for (unsigned i = 1; i <= 2; ++i) {
    final_acf(i - 1) = (Acf[i] / Acf[0]);
  }
  
  //cout << "solution = " << x << endl;
  vector<Type> Cor = RecursiveFilter(AR, nBin, final_acf);
  vector<Type> newVec(nBin);
  for (unsigned i = 0; i < nBin; ++i) {
    if (i == 0)
      newVec[i] = final_acf[0];
    else if (i == 1) {
      newVec[i] = final_acf[1];
    } else {
      newVec[i] = Cor[i];
    }
  }
  return newVec;
}
