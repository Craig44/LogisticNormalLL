/*
* The main functions returns a matrix of standardised residuals according to Francis 2014
*/
template <class Type>  
array<Type> logis_norm_stand_resids(array<Type>& obs_mat, matrix<Type>& exp_mat, vector<Type>& N,  Type sigma, vector<Type> phi, vector<Type> bin_labels, int ARMA, int centered) {
  int A = bin_labels.size();
  int Y = exp_mat.rows();
  int N_bins = exp_mat.cols();
  
  vector<Type> weights_by_year(exp_mat.rows());
  Type mean_N = N.mean();
  // Calculate a weights that are essentially scaled N
  for (int this_N = 0; this_N < N.size();++this_N) {
    weights_by_year[this_N] = sqrt(mean_N/N[this_N]);
  }
  
  // Initialise some matricies
  matrix<Type> logistic_covariance(N_bins,N_bins);
  logistic_covariance.setZero();
  matrix<Type> k_mat(N_bins - 1,N_bins);
  k_mat.setZero();
  matrix<Type> V_mat(N_bins - 1,N_bins);
  V_mat.setZero();
  logistic_covariance = covariance_logistic(sigma,phi,N_bins,ARMA);

  for (int i = 0; i < (N_bins - 1); ++i) {
    k_mat(i,i) = 1.0;
    k_mat(i, N_bins - 1) = -1.0;
  }
  matrix<Type> t_kmat = logistic_covariance * k_mat.transpose();
  
  V_mat = k_mat * t_kmat;

  if (centered == 1) {
    matrix<Type> FF(V_mat.rows(),exp_mat.cols());
    matrix<Type> H(V_mat.rows(),V_mat.cols());
    H.setOnes();
    FF.setZero();
    for (int i = 0; i < (N_bins - 1); ++i) {
      FF(i,i) = 1.0;
      FF(i, N_bins - 1) = 1.0;
      H(i,i) = 2.0;
    }
    matrix<Type> Hinv = atomic::matinv(H);
    matrix<Type> FHinv = FF.transpose() * Hinv;
    matrix<Type> Gam = FHinv * (V_mat * FHinv.transpose());
    vector<Type> diag = Gam.diagonal();
    vector<Type> sq_diag = sqrt(diag);
    matrix<Type> sdmat(exp_mat.rows(),exp_mat.cols());
    sdmat.setZero();
    for (int i = 0; i < Y; ++i)
      sdmat.row(i) = sq_diag * weights_by_year[i];
    
    array<Type> sres(exp_mat.rows(),exp_mat.cols());
    sres.setZero();
    Type Power = Type(1 / Type(obs_mat.cols()));
    for (int i = 0; i < Y; ++i) {
      Type obs_geo_mean = pow(vector<Type>(obs_mat.matrix().row(i)).prod(), Power);
      Type exp_geo_mean = pow(vector<Type>(exp_mat.row(i)).prod(),Power);
      for (int a = 0; a < N_bins; ++a) {
        sres(i,a) = (log(obs_mat(i,a) / obs_geo_mean) - log(exp_mat(i,a) / exp_geo_mean)) / sdmat(i,a);
      }
    }
    return sres;
  } else {
    vector<Type> diag = V_mat.diagonal();
    vector<Type> sq_diag = sqrt(diag);
    matrix<Type> sdmat(exp_mat.rows(),V_mat.cols());
    sdmat.setZero();
    
    for (int i = 0; i < Y; ++i)
      sdmat.row(i) = sq_diag * weights_by_year[i];
    array<Type> sres(exp_mat.rows(),V_mat.cols());
    sres.setZero();
    for (int i = 0; i < Y; ++i) {
      for (int a = 0; a < V_mat.cols(); ++a) {
        sres(i,a) = (log(obs_mat(i,a)) - obs_mat(i,(N_bins - 1))-log(exp_mat(i,a)/exp_mat(i,(N_bins - 1)))) / sdmat(i,a);
      }
    }
    return sres;
  }
  return 0;
}