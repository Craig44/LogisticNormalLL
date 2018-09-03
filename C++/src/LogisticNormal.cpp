#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <fstream>
#include <sstream>

// Include Eigen Library
#include <Eigen/Dense>

#define PI 3.14159265359
using namespace std;
using namespace Eigen;

// vector utility functions
MatrixXd readMatrix(const string& filename);
double Sum(vector<double> x);
double Mean(vector<double> x);
bool all_ones(vector<double> x);

double NLLlogistnorm(MatrixXd& obs_mat, MatrixXd& exp_mat, vector<double>& N,  double& sigma, vector<double>& phi ,bool& sepbysex,bool& sexlag, bool& robust,bool& ARMA, vector<unsigned>& bin_labels);
MatrixXd covmat_logistnorm(double& sigma, vector<double> phi,int N_bins , bool sepbysex, bool sexlag, bool ARMA, vector<unsigned> bin_labels);
MatrixXd calculate_multivariate_normal_covariance(double sigma, vector<double> phi,int N_bins, bool ARMA, bool bin_lag, vector<unsigned> bin_labels);

vector<double> ARMAacf(double AR, double MA, int nBin);
vector<double> ARMAacf(vector<double> AR, int nBin);
vector<double> GetRho(double& Phi1, int nBin);
vector<double> RecursiveFilter(vector<double> ar_coef, int nBins, vector<double> initial_vals);

double logdet(const MatrixXd& M, bool use_cholesky);
////////////////////
// Begin main function
////////////////////
int main(int argc, char* argv[]) {
  cout << "enter main" << endl;
  string task;
  int min_age;
  task = argv[1];
  min_age = atoi(argv[2]);
  bool sexed = false;
  if (task == "sex") {
    sexed = true;
  }
  bool debug = false;
  //char obs_file_name = "observed_data.txt";

  //cout << "reading matrix" << endl;

  MatrixXd obs_mat = readMatrix("observed_data.txt");
  MatrixXd exp_mat = readMatrix("expected_data.txt");
  MatrixXd error_mat = readMatrix("error_data.txt");

  //cout << "Finished reading data" << endl;


  //cout << error_mat.rows() << " cols = " << error_mat.cols() << endl;
  unsigned N_years = obs_mat.rows();
  unsigned N_bins = obs_mat.cols();


  if (obs_mat.rows() != exp_mat.rows()) {
    cerr << "observed rows != expected rows" << endl;
    return 1;
  }

  if (obs_mat.cols() != exp_mat.cols()) {
    cerr << "observed cols != expected cols" << endl;
    return 1;
  }

  vector<double> N(N_years, 0.0);
  vector<int> years;
  vector<unsigned> ages;
  int first_year = 1990;
  for (int i = 0; i < N_years; i++) {
    years.push_back(first_year + i);
    N[i] = error_mat(i,0);
  }
  unsigned N_ages = N_bins;
  if (sexed)
    N_ages = N_bins / 2;
  for (int age = 0; age < N_ages; ++age)
    ages.push_back(age + min_age);
  if (sexed) {
    for (int age = 0; age < N_ages; ++age)
      ages.push_back(age + min_age);
  }
  cout << "Nages = " << N_ages << " N bins = " << N_bins << endl;


	// Debug example values

	if (debug) {
    cout << "\n rows = " << N_years << " cols = " << N_ages << " \n";

    cout << "Observed = " << endl;
    cout << obs_mat << endl;


    cout << "Expected = " << endl;
    cout << exp_mat << endl;

    cout << "Error = " << endl;
    for (auto error : N)
      cout << error << " " ;
    cout << "\n";
	}


	bool sepbysex = false;
	bool sexlag = false;
	bool ARMA = false;
	bool robust = false;
	//double score = NLLlogistnorm(obs_mat,exp_mat,N,Sigma, phi,sepbysex, sexlag,  robust, ARMA);
	vector<double> results;
	double value;
  double sigma = 0.9;
  vector<double> phi = {0.2};
	cout << "first run " << endl;
	value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
	results.push_back(value);
	sigma = 0.2;
  cout << "Second run " << endl;
  value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
  results.push_back(value);

  phi = {0.2,0.5};
  sigma = 0.4;
  cout << "Third run " << endl;
  value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
  results.push_back(value);

  phi = {0.4,0.1};
  sigma = 0.6;
  cout << "fourth run " << endl;
  value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
  results.push_back(value);

  phi = {0.4,0.1};
  sigma = 0.6;
  ARMA = true;
  cout << "Fith run" << endl;
  value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
  results.push_back(value);

  cout << "Finished non-sexed runs" << endl;

  if (sexed) {
    // Do some additional tests for sexed data.
    sepbysex = true;
    sexlag = false;
    ARMA = true;
    phi = {0.2,0.3};
    sigma = 0.9;
    cout << "sixth run" << endl;
    value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
    results.push_back(value);

    sepbysex = false;
    sexlag = true;
    cout << "Seventh run" << endl;
    value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
    results.push_back(value);

    sexlag = false;
    value = NLLlogistnorm(obs_mat,exp_mat,N,sigma,phi,sepbysex,sexlag,robust,ARMA,ages);
    results.push_back(value);


    cout << "Finished sexed runs" << endl;

  }

  std::ofstream resultsfile("results.txt");
  resultsfile.precision(9);
  if (resultsfile.is_open())
  {
    for (auto val : results)
      resultsfile << val << " ";
  }
  resultsfile.close();

	return 0;
}

/*
 * NLL-Logistic Normal
 * @description Calculate the negative log liklihood for the logisitic-Normal distribution
 */
double NLLlogistnorm(MatrixXd& obs_mat, MatrixXd& exp_mat, vector<double>& N,  double& sigma, vector<double>& phi ,bool& sepbysex,bool& sexlag, bool& robust,bool& ARMA, vector<unsigned>& bin_labels) {
  cout << "NLLlogistnorm" << endl;
  unsigned N_bins = obs_mat.cols();
  unsigned N_years = obs_mat.rows();
  double score, mean_N;
  //int N_ages = obs_mat.cols();
  //unsigned N_years = obs_mat.rows();
  unsigned i,j;
  vector<double> weights;
  MatrixXd  Kmat(N_bins - 1,N_bins), Vmat(N_bins - 1,N_bins - 1), t_kmat, ww, V_invert, covmat;
  Kmat.setZero();
  Vmat.setZero();


  // Calculate a weights that are essentially scaled N
  score = 0.0;
  mean_N = Mean(N);
  // Calculate a weights that are essentially scaled N
  for (auto& this_N : N) {
    double temp = sqrt(mean_N/this_N);
    //cout << "weight = " << temp << endl;
    weights.push_back(temp);
  }
  // Get the covariance matrix with given parameters Sigma, phi and ARMA
  covmat = covmat_logistnorm(sigma, phi, N_bins ,sepbysex, sexlag, ARMA, bin_labels);

  std::ofstream covarfile("covariance.txt");
  if (covarfile.is_open())
  {
    covarfile << covmat << '\n';
  }
  covarfile.close();
  // Populate Kmat and Vmat
  for (i = 0; i < (N_bins - 1); ++i) {
    Kmat(i,i) = 1.0;
    Kmat(i, N_bins - 1) = -1.0;
  }

  // Calculate Vmat
  t_kmat = covmat * Kmat.transpose();
  Vmat = Kmat * t_kmat;
  //cout << "Vmat rows = " << Vmat.rows() << " cols = " << Vmat.cols() << endl;
  V_invert = Vmat.llt().solve(MatrixXd::Identity(N_bins - 1, N_bins - 1));

  std::ofstream v_invert("v_invert.txt");
  if (v_invert.is_open())
  {
    v_invert << V_invert << '\n';
  }
  v_invert.close();

  double tot_log_obs = obs_mat.array().log().sum();
  double log_det = logdet(Vmat,false);
  score = 0.5 * N_years * (N_bins - 1) * log(2 * PI) + tot_log_obs + 0.5 * N_years * log_det;

  // check if any weights deviate form 1
  if (!all_ones(weights)) {
    vector<double> log_weights;
    for (auto ww: weights)
      log_weights.push_back(log(ww));
    score += (N_bins - 1) * Sum(log_weights);
  }

  // Sweep over the the obs and create this ww object.
  ww.resize(N_years,N_bins - 1);
  for ( i = 0; i < N_years; ++i) {
    for ( j = 0; j < N_bins - 1; ++j) {
      double l_obs = obs_mat(i,j) / obs_mat(i,N_bins - 1);
      double l_exp = exp_mat(i,j) / exp_mat(i,N_bins - 1);
      ww(i, j) = log(l_obs) - log(l_exp);
    }
  }
  //cout << ww << endl;

  //cout << "\n\n\n\n" << endl;
  //cout << "----------- add robust -----------" << endl;
  // Now finish with a robustification
  if (robust) {
    // get the diaganol components of the inverted matrix
    vector<double> diag_inv;
    for (i = 0; i < V_invert.rows(); ++i)
      diag_inv.push_back(V_invert(i,i));
    for (unsigned year = 0; year < N_years; ++year) {

    }
  } else {
    double temp2;
    MatrixXd year_val, temp1;
    for (unsigned year = 0; year < N_years; ++year) {
      year_val = ww.row(year);

      temp1 = year_val * V_invert * year_val.transpose();// * year_val;

      temp2 = 0.5 / (weights[year] * weights[year]);
      //cout << temp1.size() << endl;
      score += (temp1.sum() * temp2);


    }
  }

return score;
}
/*
 * Build Up Multivariate normal covariance matrix that is structured for sexed observations it mainly calls the functions calculate_multivariate_normal_covariance();
 *
 */
MatrixXd covmat_logistnorm(double& sigma, vector<double> phi,int N_bins , bool sepbysex, bool sexlag, bool ARMA, vector<unsigned> bin_labels) {
  unsigned n_phi = phi.size();
  MatrixXd covar(N_bins,N_bins);
  covar.setZero();
  if (phi[0] == 0.0) {
    // no auto correlation
    for (unsigned diag = 0; diag < N_bins; ++ diag)
      covar(diag, diag) = sigma * sigma;
  } else if (sepbysex) {
    int N_ages = N_bins / 2; // assumes only two sexes
    for (unsigned sex = 0; sex <= 1; ++sex) {
      //cout << "nrow = " << subcov.rows() << " ncols = " << subcov.cols() << " N ages = " << N_ages << endl;
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
MatrixXd calculate_multivariate_normal_covariance(double sigma, vector<double> phi,int N_bins, bool ARMA, bool bin_lag, vector<unsigned> bin_labels) {
  MatrixXd covar(N_bins,N_bins);
  covar.setZero();
  unsigned n_phi = phi.size();

  vector<double> rho;
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
    cout << "are we here" << endl;
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

// Compute the theoretical autocorrelation function or partial autocorrelation function for an ARMA process.
// The methods used follow Brockwell & Davis (1991, section 3.3). Their equations (3.3.8) are solved for the
// autocovariances at lags 0, …, max(p, q+1), and the remaining autocorrelations are given by a recursive filter.
/**
 * @description An overloaded version of the above function, that only takes a single AR parameters.
 * @param ar_coef an AR(2) coefficient
 * @param ma_coef an MA(2) coefficien
 * @param nBin number of ages
 */
vector<double> ARMAacf(double AR, double MA, int nBin) {
  cout << "Beginning method ARMAacf() " << endl;

  // Create and declare all variables that will be used in the function
  unsigned p = 1;
  unsigned q = 1;
  vector<double> AR_coef, MA_coef,final_acf,Cor;
  AR_coef.push_back(AR);
  MA_coef.push_back(MA);

  vector<double> psi, theta, Acf;
  if (!p && !q) {
    cerr << "empty model supplied" << endl;
  }
  unsigned r = fmax(p, q + 1);

  for (unsigned i = 0; i <= (r - p); ++i) {
    AR_coef.push_back(0.0);
    p = r;
  }

  MatrixXd A(p + 1,2 * p + 1), ind(2 * p + 1,p + 1);
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


  //cout << "size of rhs " << rhs.size() << endl;
  psi.push_back(1.0);
  psi.push_back(AR + MA);
  theta.push_back(1.0);
  theta.push_back(MA);
  for (unsigned i = 0; i <= q; ++i)
    theta.push_back(0.0);

  // Declare Eigen variables
  Matrix3d A_eig; // 3f means 3x3 elements of type double
  Vector3d rhs;
  // Calculate rhs
  for (unsigned i = 0; i <= q; ++i) {
    double x1 ,x2;
    x1 = psi[0]*theta[i];
    x2 = psi[1]*theta[i + q];
    rhs(i) = Sum({x1,x2});
  }
  rhs(2) = 0.0;

  // Use the eigen library yo solve the inverse of for A with known vector B
  //vector<double> Ind;
  vector<unsigned> seq(p + 1);
  for (unsigned i = 0; i <= p; ++i) {
    seq[i] = (p - i);
  }


  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p; ++j) {
      A_eig(i,j) = A(seq[i], seq[j]);
    }
  }

  //cout << "Check A mat that we are inverting\n" << A_eig << "\n: rhs = " << rhs << endl;

  // Solve for A_eig given rhs
  Vector3d x = A_eig.colPivHouseholderQr().solve(rhs);
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
  // Print results to screen

  vector<double> newVec(nBin,0.0);
  for (unsigned i = 0; i < nBin; ++i)
    newVec[i] = Cor[i];


  return newVec;
}


/**
 * @description An overloaded version of the above function, that only takes a single AR parameters.
 * @param ar_coef an AR(2) coefficient
 * @param nBin number of ages
 */
vector<double> ARMAacf(vector<double> AR, int nBin) {
  cout << "Beginning method ARMAacf() not ARMA" << endl;
  // Create and declare all variables that will be used in the function
  unsigned p = AR.size();
  if (p > 2) {
    cerr << "This function has not been coded for more the two AR coeffs." << endl;
  }

  vector<double> rhs, psi, theta, Acf;
  unsigned r = p;
  MatrixXd A(p + 1,2 * p + 1);
  A.setZero();

  for (unsigned i = 0; i < A.rows(); ++i) {
    A(i,i) = 1.0;
    A(i,i + 1) = -AR[0];
    A(i,i + 2) = -AR[1];
  }

  rhs.assign(p + 1, 0.0);
  rhs[0] = 1.0;


  //cout << "size of rhs " << rhs.size() << endl;
  //vector<double> Ind;
  vector<unsigned> seq(p + 1);
  for (unsigned i = 0; i <= p; ++i) {
    seq[i] = (p - i);
  }

  // Declare Eigen variables
  Matrix3d A_inv; // 3f means 3x3 elements of type double
  Matrix3d A_temp; // 3f means 3x3 elements of type double
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

  //cout << "A temp " << endl;
  //cout << A_temp << endl;

  FullPivLU<Matrix3d> lu(A_temp);
  if (lu.isInvertible())
    cerr << "could not invert the matrix" << endl;

  A_inv = A_temp.inverse();

  //cout << "Check inverted matrix" << endl;
  //cout << A_inv  << endl;


  // Take the first column of the inverted matrix
  vector<double> final_acf,xx,Cor, final_Cor;
  for (unsigned i = 0; i < p + 1; ++i) {
    Acf.push_back(A_inv(i,0));
  }

  // Divide the last two elements by the first element.
  for (unsigned i = 1; i <= 2; ++i) {
    final_acf.push_back(Acf[i] / Acf[0]);
  }
  cout << "Final Acf" << endl;
  for (auto num : final_acf)
    cout << num << " ";
  cout << endl;


  // Call the recurisive filter
  Cor = RecursiveFilter(AR,nBin, final_acf);
  // Add the initial coeffiecients to the beginning of the series.
  Cor.insert(Cor.begin(),final_acf[1]);
  Cor.insert(Cor.begin(),final_acf[0]);
  // Print results to screen
  vector<double>::const_iterator first = Cor.begin();
  vector<double>::const_iterator last = Cor.begin() + nBin;
  vector<double> newVec(nBin, 0.0);
  for (unsigned i = 0; i < nBin; ++i)
    newVec[i] = Cor[i];

  cout << "\n\n\n\n";
  cout << newVec.size() << " vs " << Cor.size() << endl;
  cout << "\n\n\n\n";

  for (auto num : newVec)
    cout << num << " ";
  cout << endl;
  for (auto num : Cor)
    cout << num << " ";
  cout << endl;
  return newVec;
}

/**
 * This method is called at in the CalculateCovarianceLogisiticNormal() method to calculate the auto-regressive vector Rho
 * @param Phi1 an AR(1) coefficient
 * @param nBin number of ages
 */
vector<double> GetRho(double& Phi1, int nBin) {
  cout <<  "entering GetRho()" << endl;
  //calculation of AR(1) acf for  LN2
  vector<double> rho(nBin - 1, 0.0);

  for(int i = 1; i <= nBin - 1; i++) {
    rho[i - 1]= pow(Phi1,(double)i);
    //cout << rho[i - 1] << " ";
  }
  cout <<  "\n\n\n";

  return rho;
}


/**
 * @param ar_coef an AR(2) coefficient
 * @param nBin number of ages
 * @param initial_vals initial vals
 */
vector<double>  RecursiveFilter(vector<double> ar_coef, int nBins, vector<double> initial_vals) {
  vector<double> store_vec(nBins + 1,0.0);
  if (ar_coef.size() > 2) {
    cerr <<  "RecursiveFilter(): has not been coded for more than 2 AR coeffiecients" << endl;
  }

  store_vec[0] = initial_vals[1];
  store_vec[1] = initial_vals[0];
  vector<double> result(store_vec.size() - 1);

  for (unsigned i = 1; i < nBins + 1; ++i) {
    if (i == 1) {
      store_vec[i] =   store_vec[i - 1] *ar_coef[0]  + store_vec[i] *  ar_coef[1];
    } else {
      store_vec[i] = store_vec[i - 1] *  ar_coef[0] + store_vec[i - 2] * ar_coef[1];
    }
    //cout << "value = " << store_vec[i];
  }
  // remove the first value
  for (unsigned i = 1; i < store_vec.size(); ++i)
    result[i - 1] = store_vec[i];

  return result;
}


//////////////////////
// Utility Functions
//////////////////////
/**
 * @param Character for a filename to read
 * This function reads in a matrix and returns an eigen matrix
 */

MatrixXd readMatrix(const string& filename) {
  ifstream infile;
  infile.open(filename);
  unsigned nrows = 0;
  unsigned ncols = 1;
  if (infile.is_open()) {
      //cout << "start reading the file" << endl;
      string current_line;

      while(getline(infile, current_line)) {
        if (nrows == 0) {
          // count cols
          size_t pos = 0;
          string token;
          string delimiter = " ";
          while ((pos = current_line.find(delimiter)) != std::string::npos) {
            token = current_line.substr(0, pos);
            //std::cout << token << std::endl;
            current_line.erase(0, pos + delimiter.length());
            ++ncols;
          }
        }
        ++nrows;
      }
      infile.close();
  } else {
      cout << "could not open the file" << endl;
  }
  //cout << "matrix had n-rows = " << nrows << " and cols = " << ncols << endl;

  /// Now create a matrix and save the informations
  MatrixXd return_mat(nrows,ncols);
  infile.open(filename);
  string current_line;
  unsigned row = 0;
  unsigned col = 0;

  while(getline(infile, current_line)) {
      stringstream stream(current_line);
      col = 0;
      while(! stream.eof()) {
      stream >> return_mat(row,col);
      ++col;
      }
      ++row;
  }
  //cout << "mat" << endl;
  return return_mat;
}

double logdet(const MatrixXd& M, bool use_cholesky) {
  double ld = 0;
  if (use_cholesky) {
    LLT<Matrix<double,Dynamic,Dynamic>> chol(M);
    auto& U = chol.matrixL();
    for (unsigned i = 0; i < M.rows(); ++i)
      ld += log(U(i,i));
    ld *= 2;
  } else {
    PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu(M);
    auto& LU = lu.matrixLU();
    double c = lu.permutationP().determinant(); // -1 or 1
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const auto& lii = LU(i,i);
      if (lii < double(0)) c *= -1;
      ld += log(abs(lii));
    }
    ld += log(c);
  }
  return ld;
}
// Utilty funciton for summing
double Sum(vector<double> x) {
  double Total = 0.0;
    for (unsigned i = 0; i< x.size(); ++i) {
     Total += x[i];
    }
    return Total;
}

double Mean(vector<double> x) {
  double Total = 0.0;
    for (unsigned i = 0; i< x.size(); ++i) {
     Total += x[i];
    }
    return Total / (double)x.size();
}

bool all_ones(vector<double> x) {
  for(auto num : x) {
    if (num != 1.0)
      return false;
  }
  return true;
}
