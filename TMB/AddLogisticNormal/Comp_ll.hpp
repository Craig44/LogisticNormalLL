/*
 * NLL-Logistic Normal
 * @description Calculate the negative log liklihood for the logisitic-Normal distribution
 */
template <class Type>  
Type NLLlogistnorm(matrix<Type> obs_mat, matrix<Type> exp_mat, vector<Type> N,  Type sigma, vector<Type> phi ,bool& sepbysex,bool& sexlag, bool& robust,bool& ARMA, vector<Type> bin_labels) {
  unsigned N_bins = obs_mat.cols();
  unsigned N_years = obs_mat.rows();
  Type score, mean_N;
  //int N_ages = obs_mat.cols();
  //unsigned N_years = obs_mat.rows();
  unsigned i;
  vector<Type> weights(N.size());
  matrix<Type>  Kmat(N_bins - 1,N_bins), Vmat(N_bins - 1,N_bins - 1), t_kmat, ww, V_invert, covmat;
  Kmat.setZero();
  Vmat.setZero();


  // Calculate a weights that are essentially scaled N
  score = 0.0;
  mean_N = N.mean();
  // Calculate a weights that are essentially scaled N
  for (int this_N = 0; this_N < N.size();++this_N) {
    Type temp = sqrt(mean_N/N[this_N]);
    //cout << "weight = " << temp << endl;
    weights[this_N] = temp;
  }
  
  // Get the covariance matrix with given parameters Sigma, phi and ARMA
  covmat = covmat_logistnorm(sigma, phi, N_bins ,sepbysex, sexlag, ARMA, bin_labels);
  /*
  // Populate Kmat and Vmat
  for (i = 0; i < (N_bins - 1); ++i) {
    Kmat(i,i) = 1.0;
    Kmat(i, N_bins - 1) = -1.0;
  }

  // Calculate Vmat
  t_kmat = covmat * Kmat.transpose();
  Vmat = Kmat * t_kmat;
  V_invert = Vmat.inverse();

  Type tot_log_obs = obs_mat.array().log().sum();
  Type log_det = log(Vmat.determinant());
  score = 0.5 * N_years * (N_bins - 1) * log(2 * PI) + tot_log_obs + 0.5 * N_years * log_det;
  
  // check if any weights deviate form 1
  if (!all_ones(weights)) {
    vector<Type> log_weights;
    for (auto ww: weights)
      log_weights.push_back(log(ww));
    score += (N_bins - 1) * Sum(log_weights);
  }

  // Sweep over the the obs and create this ww object.
  ww.resize(N_years,N_bins - 1);
  for ( i = 0; i < N_years; ++i) {
    for ( j = 0; j < N_bins - 1; ++j) {
      Type l_obs = obs_mat(i,j) / obs_mat(i,N_bins - 1);
      Type l_exp = exp_mat(i,j) / exp_mat(i,N_bins - 1);
      ww(i, j) = log(l_obs) - log(l_exp);
    }
  }
  //cout << ww << endl;

  
  // Now finish with a robustification
  if (robust) {
    // get the diaganol components of the inverted matrix
    vector<Type> diag_inv;
    for (i = 0; i < V_invert.rows(); ++i)
      diag_inv.push_back(V_invert(i,i));
    for (unsigned year = 0; year < N_years; ++year) {

    }
  } else {
    Type temp2;
    MatrixXd year_val, temp1;
    for (unsigned year = 0; year < N_years; ++year) {
      year_val = ww.row(year);

      temp1 = year_val * V_invert * year_val.transpose();// * year_val;

      temp2 = 0.5 / (weights[year] * weights[year]);
      //cout << temp1.size() << endl;
      score += (temp1.sum() * temp2);
    }
  }
   */

  return 0.0;
}

/*
 * Build Up Multivariate normal covariance matrix that is structured for sexed observations it mainly calls the functions calculate_multivariate_normal_covariance();
 *
 */
template <class Type> 
matrix<Type> covmat_logistnorm(Type sigma, vector<Type> phi,unsigned N_bins , bool sepbysex, bool sexlag, bool ARMA, vector<Type> bin_labels) {
  unsigned n_phi = phi.size();
  matrix<Type> covar(N_bins,N_bins);
  covar.setZero();
  if (phi[0] == 0.0) {
    // no auto correlation
    for (unsigned diag = 0; diag < N_bins; ++ diag)
      covar(diag, diag) = sigma * sigma;
  } else if (sepbysex) {
    int N_ages = N_bins / 2; // assumes only two sexes
    for (unsigned sex = 0; sex <= 1; ++sex) {
      covar.block(sex * N_ages, sex * N_ages, N_ages, N_ages) = calculate_multivariate_normal_covariance(sigma, phi, N_ages, ARMA, false, bin_labels);
    }
  } else if (sexlag) {
    covar = calculate_multivariate_normal_covariance(sigma, phi, N_bins, ARMA, true, bin_labels);
  } else {
    covar = calculate_multivariate_normal_covariance(sigma, phi, N_bins, ARMA, false, bin_labels);
  }
  return covar;
}


/*
 * @function calculate_multivariate_normal_covariance
 * @description a generic function that calcualtes and arbituary multivariate covariance function with some correltion
 * @return a covariance dictated by the parameters.
 *
 */
template <class Type> 
matrix<Type> calculate_multivariate_normal_covariance(Type sigma, vector<Type> phi,int N_bins, bool ARMA, bool bin_lag, vector<Type> bin_labels) {
  matrix<Type> covar(N_bins,N_bins);
  covar.setZero();
  unsigned n_phi = phi.size();

  vector<Type> rho;
  if (bin_lag) {
    if (n_phi == 1) {
      rho = GetRho(phi[0],N_bins + 1);
    } else if (n_phi == 2 && ARMA) {
      rho = ARMAacf(phi[0],phi[1],N_bins + 1);
    } else {
      rho = ARMAacf(phi,N_bins + 1);
    }

    for (int diag = 0; diag < bin_labels.size(); ++ diag) {
      covar(diag,diag) = 1.0 * sigma * sigma;
    }
    for (int diag = 0; diag < N_bins; ++ diag) {
      for (int row = 0; row < bin_labels.size(); ++ row) {
        for (int col = 0; col < bin_labels.size(); ++ col) {
          if ((fabs(bin_labels[row] - bin_labels[col]) + 1) == (diag + 1)) {
            covar(row, col) = rho[diag] * sigma * sigma;
            covar(col, row) = rho[diag] * sigma * sigma; // make it symetric
          }
          if (row == col)
            covar(col, row) = sigma * sigma;
        }
      }
    }
  } else {
    if (n_phi == 1) {
      rho = GetRho(phi[0],N_bins);
    } else if (n_phi == 2 && ARMA) {
      rho = ARMAacf(phi[0],phi[1],N_bins);
    } else {
      rho = ARMAacf(phi,N_bins);
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
 * @description Calcualte the logdeterminant of a SPD matrix
 * @param M the matrix to perform the operation on
 * @param use_cholesky bool whether to use the cholesky decomposition
 * @return a singl value of the log determinant of the matrix
 */
/*
template <class Type> 
Type logdet(const matrix<Type>& M, bool use_cholesky) {
  Type ld = 0;
  if (use_cholesky) {
    LLT<matrix<Type,Dynamic,Dynamic>> chol(M);
    auto& U = chol.matrixL();
    for (unsigned i = 0; i < M.rows(); ++i)
      ld += log(U(i,i));
    ld *= 2;
  } else {
    PartialPivLU<matrix<Type,Dynamic,Dynamic>> lu(M);
    auto& LU = lu.matrixLU();
    Type c = lu.permutationP().determinant(); // -1 or 1
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const auto& lii = LU(i,i);
      if (lii < Type(0)) c *= -1;
      ld += log(abs(lii));
    }
    ld += log(c);
  }
  return ld;
}
 */
/*
 * --------------------------------------
 *        Auto-Correlation functions
 * --------------------------------------
*/

/**
 * @description Calcualte the ARMA correlation structure in the matrix
 * @param ar_coef an AR(2) coefficient
 * @param ma_coef an MA(2) coefficien
 * @param nBin number of ages
 */
template <class Type> 
vector<Type> ARMAacf(Type AR, Type MA, int nBin) {
  // Create and declare all variables that will be used in the function
  unsigned p = 1;
  unsigned q = 1;
  vector<Type> AR_coef, MA_coef,final_acf,Cor;
  AR_coef.push_back(AR);
  MA_coef.push_back(MA);

  vector<Type> psi, theta, Acf;
  if (!p && !q) {
    std::cerr << "empty model supplied" << std::endl;
  }
  unsigned r = fmax(p, q + 1);

  for (unsigned i = 0; i <= (r - p); ++i) {
    AR_coef.push_back(0.0);
    p = r;
  }

  matrix<Type> A(p + 1,2 * p + 1), ind(2 * p + 1,p + 1);
  A.setZero();
  ind.setZero();
  for (unsigned i = 0; i < ind.rows(); ++i) {
    for (unsigned j = 0; j < ind.cols(); ++j) {
      ind(i,0) = i + 1;
      ind(i, 1) = (i + 1) + (j + 1) - (i + 1);
    }
  }
  for (unsigned i = 0; i < A.rows(); ++i) {
    A(i ,i) = 1.0;
    A(i, i + 1) = -AR_coef[0];
  }

  A(0 ,A.rows() - 1) = -AR_coef[1];
  A(A.rows() - 1, 0) = -AR_coef[1];
  A(A.rows() - 1, 1) = -AR_coef[0];


  psi.push_back(1.0);
  psi.push_back(AR + MA);
  theta.push_back(1.0);
  theta.push_back(MA);
  for (unsigned i = 0; i <= q; ++i)
    theta.push_back(0.0);

  // Declare Eigen variables
  matrix<Type> A_eig(3,3); // 3f means 3x3 elements of type Type
  vector<Type> rhs(3,0.0);
  // Calculate rhs
  for (unsigned i = 0; i <= q; ++i) {
    Type x1 ,x2;
    x1 = psi[0]*theta[i];
    x2 = psi[1]*theta[i + q];
    rhs(i) = x1 + x2;
  }
  rhs(2) = 0.0;

  // Use the eigen library yo solve the inverse of for A with known vector B
  //vector<Type> Ind;
  vector<unsigned> seq(p + 1);
  for (unsigned i = 0; i <= p; ++i) {
    seq(i) = (p - i);
  }


  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p; ++j) {
      A_eig(i,j) = A(seq[i], seq[j]);
    }
  }


  // Solve for A_eig given rhs
  matrix<Type> x = A_eig.colPivHouseholderQr().solve(rhs);
  //cout << "solution = " << x << endl;
  
  
  
  // Divide the last two elements by the first element.
  //cout << "Find the crash" << endl;

  for (unsigned i = 1; i <= 2; ++i) {
    final_acf.push_back(x(i) / x(0));
  }


  Cor = RecursiveFilter(AR_coef, nBin, final_acf);

  // Add the initial coeffiecients to the beginning of the series.
  Cor.insert(Cor.begin(), final_acf[1]);
  Cor.insert(Cor.begin(), final_acf[0]);
   vector<Type> newVec(nBin,0.0);
   for (unsigned i = 0; i < nBin; ++i)
  newVec[i] = Cor[i];
  
  
  return newVec;
}
/**
 * @description An overloaded version of the above function, that only takes a single AR parameters.
 * @param ar_coef an AR(2) coefficient
 * @param nBin number of ages
 */
template <class Type> 
vector<Type> ARMAacf(vector<Type> AR, int nBin) {
  // Create and declare all variables that will be used in the function
  unsigned p = AR.size();
  if (p > 2) {
    std::cerr << "This function has not been coded for more the two AR coeffs.";
  }

  vector<Type> rhs, psi, theta, Acf;
  unsigned r = p;
  matrix<Type> A(p + 1,2 * p + 1);
  A.setZero();

  for (unsigned i = 0; i < A.rows(); ++i) {
    A(i,i) = 1.0;
    A(i,i + 1) = -AR[0];
    A(i,i + 2) = -AR[1];
  }

  rhs.assign(p + 1, 0.0);
  rhs[0] = 1.0;

  vector<unsigned> seq(p + 1);
  for (unsigned i = 0; i <= p; ++i) {
    seq(i) = (p - i);
  }
  // Declare Eigen variables
  matrix<Type> A_inv(3,3); // 3f means 3x3 elements of type Type
  matrix<Type> A_temp(3,3); // 3f means 3x3 elements of type Type
  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p ; ++j) {
      if (j == 2)
        A_temp(i,j) = A(i,j);
      else
        A_temp(i,j) = A(i,j) + A(i,2 * p  - j);
    }
  }

  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p ; ++j) {
      A_temp(i,j) =  A_temp(seq[i],seq[j]);
    }
  }
  // the bodge
  A_temp(1,2) = 0.0;

  A_inv = A_temp.inverse();

  // Take the first column of the inverted matrix
  vector<Type> final_acf,xx,Cor, final_Cor;
  for (unsigned i = 0; i < p + 1; ++i) {
    Acf.push_back(A_inv(i,0));
  }

  // Divide the last two elements by the first element.
  for (unsigned i = 1; i <= 2; ++i) {
    final_acf.push_back(Acf[i] / Acf[0]);
  }


  // Call the recurisive filter
  xx.assign(nBin - p, 0.0);
  Cor = RecursiveFilter(AR,nBin, final_acf);
  // Add the initial coeffiecients to the beginning of the series.
  Cor.insert(Cor.begin(),final_acf[1]);
  Cor.insert(Cor.begin(),final_acf[0]);
  // Print results to screen
  vector<Type> newVec(nBin,0.0);
  for (unsigned i = 0; i < nBin; ++i)
    newVec[i] = Cor[i];
  
  return newVec;
}
/**
 * This method is called at in the CalculateCovarianceLogisiticNormal() method to calculate the auto-regressive vector Rho
 * @param Phi1 an AR(1) coefficient
 * @param nBin number of ages
 */
template <class Type> 
vector<Type> GetRho(Type& Phi1, int nBin) {
  //calculation of AR(1) acf for  LN2
  vector<Type> rho(nBin - 1, 0.0);

  for(int i = 1; i <= nBin - 1; i++) {
    rho[i - 1]= pow(Phi1,(Type)i);
    //cout << rho[i - 1] << " ";
  }

  return rho;
}

/**
 * @RecursiveFilter 
 * @param ar_coef an AR(2) coefficient
 * @param nBin number of ages
 * @param initial_vals initial vals
 */
template <class Type> 
vector<Type>  RecursiveFilter(vector<Type> ar_coef, int nBins, vector<Type> initial_vals) {
  vector<Type> store_vec(nBins + 1,0.0);
  
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

template <class Type> 
bool all_ones(vector<Type> x) {
  for(unsigned i = 0; i < x.size(); ++i) {
    if (x[i] != 1.0)
      return false;
  }
  return true;
}

