#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <fstream>

// Tapping into Casal2's boost library
#include <boost\numeric\ublas\matrix.hpp>
#include <boost\numeric\ublas\lu.hpp>
#include <boost\numeric\ublas\io.hpp>
#include <boost\numeric\ublas\storage.hpp>
#include <Eigen/Dense>

#define PI 3.14159265359
using namespace std;
//using namespace ublas = boost::numeric::ublas;

// Evaluate the likelihood contribution for the logistic normal, with the following parameters.
/*
 *@param obs Nbin-vector, or a Nyear x Nbin matrix of Nyear observed Nbin-vectors of proportions
 *@param exp Nbin-vector, or a Nyear x Nbin matrix of Nyear expected Nbin-vectors of proportions
 *@param N a Nyear vector of effective sample sizes
 *@param sigma an estimable standard deviation parameter
 *@param phi  auto correlation parameter
 *@param covmat a user defined covariance matrix with the same dimensions as obs
 *@param sepbysex ignore the rest of the bool switches they are for dealing with sexed populations
 *
 *@return a negative log likelihood score for exp given observations obs and parameters parameters, sigma and phi
 */
double NLLlogistnorm(vector<vector<double>> obs,vector<vector<double>> exp ,vector<double> N, double sigma, vector<double> phi ,bool sepbysex,bool sexlag, bool robust,bool ARMA);

// Generate a set of simulated observations that follow a logistic-normal distribution that
vector<double> rlogistnorm(vector<double> expprop, double sigma, vector<double> phi, vector<vector<double>> covmat, bool ARMA);


// @param sigmaeither a single number or a vector of the same length as binnam, defining the s.d. of X
// @param phi a 1- or 2-vector defining the correlation structure of X
// @param sepbysex sexlag - options for sexed compositions (see function NLLlogistnorm for details)
// @param binnam vector of strings that contain the bin labels
// @param ARMA if T and if phi is of length 2 than X is assumed to have the correlation matrix of an ARMA(1,1) process, where phi[1] and phi[2] are the AR and MA parameters, respectively
vector<vector<double> > covmat_logistnorm(double sigma, vector<double> phi,unsigned N_ages , bool sepbysex, bool sexlag, bool ARMA);  // Forward declaration


// @param AR_coef numeric vector of AR coefficients
// @param MA_coef numeric vector of MA coefficients
// @param lag_max Maximum lag required. Defaults to max(p, q+1), where p, q are the numbers of AR and MA terms respectively.
// @param pacf Should the partial autocorrelations be returned?
// Thwo version of this funciton, one with a two AR and another with a AR and MA coeffecient
// For information on this algorithm, check the ARMAacf() function in R that was used as the basis of testing the algorithm.
vector<double> ARMAacf(double AR, double MA, int nBin);
vector<double> ARMAacf(vector<double> AR, int nBin);
vector<double> GetRho(double& Phi1, int nBin);
vector<double>  RecursiveFilter(vector<double> ar_coef, int nBins, vector<double> initial_vals);

// Forward declarations for utility funcitons
double Sum(vector<double> x);
double Mean(vector<double> x);
bool all_ones(vector<double> x);
vector<vector<double> > mat_multiply(vector<vector<double> >& x,vector<vector<double> >& y);
vector<vector<double> > mat_multiply(vector<double> & x,vector<vector<double> >& y);
vector<vector<double> > mat_multiply(vector<vector<double> > & x,vector<double>& y);

// transpose matrix x
vector<vector<double> > t(vector<vector<double> >& x);
vector<vector<double>> log_mat(vector<vector<double> >& x);
double Sum_mat(vector<vector<double> >& x);

vector<double> elem_prod(vector<double> x,vector<double> y);
bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);
double det_fast(const boost::numeric::ublas::matrix<double>& matrix);
////////////////////
// Begin main function
////////////////////
int main() {
	// Debug example values
	double Sigma = 0.4;
	vector<double> phi = {0.2};
	vector<string> Binnam = {"1", "2", "3", "4", "5","6","7","8", "9", "10", "11", "12","15","16","17","18", "19",};
	vector<vector<double> > obs_mat = {{0.0019, 0.0215, 0.0275, 0.0158, 0.0568, 0.0499, 0.1001, 0.2629, 0.1331, 0.1034, 0.1294, 0.0231, 0.0272, 0.0002, 0.0132, 0.0110, 0.0115, 0.0002, 0.0115},
	    {0.0213, 0.0350, 0.0420, 0.0101, 0.0255, 0.0449, 0.0950, 0.1354, 0.0619, 0.2855, 0.0864, 0.0746, 0.0508, 0.0180, 0.0129, 0.0002, 0.0002, 0.0002, 0.0002},
	    {0.0104, 0.0435, 0.0641, 0.0625, 0.0576, 0.0398, 0.0943, 0.0586, 0.0479, 0.0858, 0.1897, 0.0945, 0.0736, 0.0319, 0.0157, 0.0002, 0.0063, 0.0064, 0.0171},
	    {0.0310, 0.1306, 0.1041, 0.0186, 0.0106, 0.0459, 0.0345, 0.0644, 0.0518, 0.0281, 0.0358, 0.1525, 0.0992, 0.0670, 0.0502, 0.0397, 0.0255, 0.0052, 0.0052},
	    {0.0785, 0.1866, 0.1408, 0.1011, 0.0845, 0.0440, 0.0274, 0.0741, 0.0488, 0.0475, 0.0238, 0.0115, 0.0170, 0.0248, 0.0034, 0.0130, 0.0099, 0.0015, 0.0617},
	    {0.0451, 0.2763, 0.1077, 0.0613, 0.0739, 0.0454, 0.0698, 0.0293, 0.0488, 0.0517, 0.0454, 0.0305, 0.0260, 0.0183, 0.0040, 0.0221, 0.0137, 0.0002, 0.0304},
	    {0.0172, 0.2318, 0.1036, 0.0976, 0.0880, 0.0648, 0.0678, 0.0740, 0.0429, 0.0434, 0.0316, 0.0309, 0.0321, 0.0072, 0.0173, 0.0036, 0.0055, 0.0165, 0.0242},
	    {0.0319, 0.1448, 0.1475, 0.1015, 0.1384, 0.1486, 0.0737, 0.0343, 0.0426, 0.0312, 0.0414, 0.0153, 0.0060, 0.0107, 0.0002, 0.0002, 0.0093, 0.0080, 0.0147},
	    {0.0083, 0.0885, 0.0926, 0.0979, 0.1544, 0.1528, 0.0725, 0.1093, 0.0293, 0.0224, 0.0347, 0.0257, 0.0281, 0.0161, 0.0149, 0.0086, 0.0078, 0.0002, 0.0360},
	    {0.0263, 0.1047, 0.1482, 0.1406, 0.0966, 0.1083, 0.0926, 0.0896, 0.0230, 0.0781, 0.0389, 0.0104, 0.0159, 0.0047, 0.0020, 0.0023, 0.0010, 0.0010, 0.0157},
	    {0.1060, 0.1153, 0.0800, 0.0675, 0.0758, 0.1056, 0.0984, 0.1125, 0.0508, 0.0681, 0.0557, 0.0131, 0.0246, 0.0002, 0.0091, 0.0096, 0.0059, 0.0002, 0.0016},
	    {0.0993, 0.1142, 0.1035, 0.0734, 0.0911, 0.0912, 0.0900, 0.0774, 0.0545, 0.0644, 0.0352, 0.0429, 0.0205, 0.0126, 0.0160, 0.0026, 0.0059, 0.0002, 0.0050},
	    {0.1408, 0.2712, 0.1572, 0.0636, 0.0537, 0.0767, 0.0348, 0.0312, 0.0243, 0.0243, 0.0162, 0.0367, 0.0127, 0.0085, 0.0157, 0.0091, 0.0002, 0.0013, 0.0220},
	    {0.2373, 0.3437, 0.1528, 0.0439, 0.0199, 0.0287, 0.0344, 0.0199, 0.0266, 0.0180, 0.0111, 0.0216, 0.0077, 0.0234, 0.0002, 0.0035, 0.0002, 0.0002, 0.0070},
	    {0.0139, 0.1752, 0.1619, 0.1173, 0.0827, 0.0573, 0.0735, 0.0820, 0.0387, 0.0323, 0.0351, 0.0228, 0.0177, 0.0141, 0.0263, 0.0188, 0.0045, 0.0049, 0.0208},
	    {0.0347, 0.2180, 0.1783, 0.0999, 0.0856, 0.0776, 0.0659, 0.0656, 0.0432, 0.0180, 0.0177, 0.0203, 0.0173, 0.0053, 0.0043, 0.0179, 0.0095, 0.0129, 0.0081}};

  vector<vector<double> > exp_mat = {
{0.03014, 0.04666, 0.04960, 0.05579, 0.07213, 0.05041, 0.10041, 0.17720, 0.12866, 0.08687, 0.06243, 0.03883, 0.02502, 0.01627, 0.01287, 0.01139, 0.00950, 0.00740, 0.01840},
{0.04074, 0.07713, 0.06410, 0.04238, 0.04247, 0.05165, 0.06134, 0.05454, 0.09144, 0.13744, 0.11152, 0.07336, 0.05044, 0.03180, 0.01967, 0.01266, 0.00918, 0.00751, 0.02063},
{0.05752, 0.09241, 0.07116, 0.05653, 0.03997, 0.03988, 0.04823, 0.05518, 0.05298, 0.08222, 0.11611, 0.09747, 0.06459, 0.04320, 0.02732, 0.01664, 0.01055, 0.00737, 0.02067},
{0.06083, 0.12746, 0.08396, 0.06245, 0.05016, 0.03722, 0.03653, 0.04322, 0.04814, 0.04879, 0.07072, 0.09453, 0.08128, 0.05458, 0.03561, 0.02246, 0.01351, 0.00840, 0.02015},
{0.07894, 0.15535, 0.12990, 0.11899, 0.09094, 0.07396, 0.06944, 0.05721, 0.04365, 0.03168, 0.02272, 0.01694, 0.01306, 0.01120, 0.01081, 0.01089, 0.01116, 0.01183, 0.04132},
{0.06990, 0.16744, 0.13393, 0.10943, 0.10010, 0.07856, 0.06359, 0.05771, 0.04743, 0.03581, 0.02585, 0.01836, 0.01346, 0.01027, 0.00862, 0.00807, 0.00795, 0.00796, 0.03555},
{0.05178, 0.15175, 0.14689, 0.11589, 0.09694, 0.08782, 0.06960, 0.05570, 0.04910, 0.04015, 0.03008, 0.02155, 0.01515, 0.01093, 0.00822, 0.00674, 0.00612, 0.00587, 0.02971},
{0.05546, 0.11506, 0.13669, 0.12861, 0.10482, 0.08821, 0.07807, 0.06189, 0.04899, 0.04208, 0.03415, 0.02541, 0.01805, 0.01256, 0.00892, 0.00660, 0.00528, 0.00464, 0.02452},
{0.05716, 0.12394, 0.10599, 0.11833, 0.11370, 0.09422, 0.07845, 0.06761, 0.05352, 0.04200, 0.03521, 0.02827, 0.02091, 0.01471, 0.01013, 0.00707, 0.00513, 0.00400, 0.01964},
{0.07900, 0.12371, 0.10817, 0.09195, 0.10159, 0.09845, 0.08183, 0.06718, 0.05666, 0.04476, 0.03481, 0.02850, 0.02260, 0.01660, 0.01157, 0.00787, 0.00540, 0.00384, 0.01550},
{0.10140, 0.16201, 0.10433, 0.08622, 0.07639, 0.08293, 0.07992, 0.06620, 0.05368, 0.04444, 0.03493, 0.02692, 0.02154, 0.01683, 0.01227, 0.00846, 0.00568, 0.00383, 0.01203},
{0.10106, 0.20025, 0.13008, 0.08257, 0.06809, 0.06184, 0.06533, 0.06227, 0.05146, 0.04130, 0.03358, 0.02621, 0.01999, 0.01565, 0.01203, 0.00869, 0.00593, 0.00393, 0.00975},
{0.09751, 0.20093, 0.16020, 0.10266, 0.06709, 0.05474, 0.04976, 0.05085, 0.04788, 0.03948, 0.03138, 0.02507, 0.01939, 0.01463, 0.01121, 0.00847, 0.00606, 0.00409, 0.00861},
{0.08993, 0.19600, 0.16367, 0.12598, 0.08360, 0.05555, 0.04439, 0.03998, 0.03967, 0.03688, 0.03033, 0.02386, 0.01875, 0.01433, 0.01069, 0.00802, 0.00595, 0.00420, 0.00822},
{0.07011, 0.15630, 0.15514, 0.13439, 0.11251, 0.08908, 0.06132, 0.04084, 0.03132, 0.02747, 0.02589, 0.02336, 0.01904, 0.01467, 0.01114, 0.00828, 0.00600, 0.00432, 0.00882},
{0.03833, 0.15238, 0.13923, 0.13306, 0.11866, 0.10056, 0.07909, 0.05446, 0.03615, 0.02724, 0.02345, 0.02156, 0.01914, 0.01549, 0.01181, 0.00881, 0.00645, 0.00460, 0.00951}};
	vector<double> N = {19, 21, 30, 36, 58, 46, 52, 38, 30, 40, 51, 49, 59, 45, 49, 60};
	bool sepbysex = false;
	bool sexlag = false;
	bool ARMA = false;
	bool robust = false;
	double score = NLLlogistnorm(obs_mat,exp_mat,N,Sigma, phi,sepbysex, sexlag,  robust, ARMA);



	vector<vector<double>> covmat = covmat_logistnorm(Sigma, phi,19,sepbysex, sexlag, ARMA);

	// End function
	return 0;
}

// Define the Negative loglikeihood score for a logistic normal distributed comp data.
double NLLlogistnorm(vector<vector<double>> obs,vector<vector<double>> exp,vector<double> N,  double sigma, vector<double> phi ,bool sepbysex,bool sexlag, bool robust,bool ARMA) {
  double score, mean_N, N_years, N_ages;
  unsigned i,j,k;
  vector<double> weights;
  vector<vector<double>> covmat, Kmat, Vmat,t_kmat,log_obs,ww,V_invert;
  N_years = obs.size();
  N_ages = obs[0].size();
  // initialise Kmat and Vmat;
  Kmat.resize(N_ages, vector<double>(N_ages,0.0));
  Vmat.resize(N_ages, vector<double>(N_ages,0.0));
  score = 0.0;
  mean_N = Mean(N);
  cout << "mean N = " << mean_N << ", no years = " << N_years << ", no ages = " << N_ages << endl;
  // Calculate a weights that are essentially scaled N
  for (auto this_N : N) {
    double temp = sqrt(mean_N/this_N);
    cout << "weight = " << temp << endl;
    weights.push_back(temp);
  }
  // Get the covariance matrix with given parameters Sigma, phi and ARMA
  covmat = covmat_logistnorm(sigma, phi, N_ages,sepbysex, sexlag, ARMA);

  cout << "printing covar = " << endl;
  for (k = 0; k < covmat.size(); ++k) {
    for (j = 0; j < covmat[k].size(); ++j)
      cout << covmat[k][j] << " ";
    cout << endl;
  }

  // Populate Kmat and Vmat
  for (i = 0; i < (N_ages - 1); ++i) {
    Kmat[i][i] = 1.0;
    Kmat[i][N_ages - 1] = -1.0;
  }

  // Calculate Vmat
  t_kmat = t(Kmat);
  t_kmat= mat_multiply(covmat,t_kmat);

  Vmat = mat_multiply(Kmat, t_kmat);

  cout << "Size of Kmat =" <<Kmat.size() << " size of t_kmat= " << t_kmat.size() <<  endl;

  cout << "printing Vmat = " << endl;
  boost::numeric::ublas::matrix<double> ublasV_mat(N_ages - 1, N_ages - 1);
  boost::numeric::ublas::matrix<double> inv_Vmat(N_ages - 1, N_ages - 1);
  for (k = 0; k < Vmat.size() - 1; ++k) {
    for (j = 0; j < Vmat[k].size() - 1; ++j) {
      cout << Vmat[k][j] << " ";
      ublasV_mat(k,j) = Vmat[k][j];
    }
    cout << endl;
  }
  bool inverted = InvertMatrix(ublasV_mat,inv_Vmat);
  if (!inverted) {
    cerr << "could not invert matrix, exiting program";
    return 0;
  }
  // Need to convert inverse back into a function that will be handled by all my shitty functions... how fustrating, this is where templating will be useful
  V_invert.resize(inv_Vmat.size1(), vector<double>(inv_Vmat.size2()));
  cout << inv_Vmat.size1() << " " << inv_Vmat.size2() << endl;
  for (i = 0; i < inv_Vmat.size1(); ++i) {
    for (j = 0; j < inv_Vmat.size2(); ++j) {
      V_invert[i][j] = inv_Vmat(i,j);
    }
    cout << endl;
  }

  log_obs = log_mat(obs);

  double tot_log_obs = Sum_mat(log_obs);
  score = 0.5 * N_years * (N_ages - 1) * log(2 * PI) + tot_log_obs + 0.5 * N_years * log(det_fast(ublasV_mat));
  cout << "total log obs = " << tot_log_obs << " determinant of inverted mat = " << det_fast(ublasV_mat) << " initial score = " <<  score << endl;

  // check if any weights deviate form 1
  if (!all_ones(weights)) {
    vector<double> log_weights;
    for (auto ww: weights)
      log_weights.push_back(log(ww));
    score += (N_ages - 1) * Sum(log_weights);
  }

  // Sweep over the the obs and create this ww object.
  ww.resize(N_years, vector<double>(N_ages - 1));
  for ( i = 0; i < N_years; ++i) {
    for ( j = 0; j < N_ages - 1; ++j) {
      double l_obs = obs[i][j] / obs[i][N_ages - 1];
      double l_exp = exp[i][j] / exp[i][N_ages - 1];
      ww[i][j] = log(l_obs) - log(l_exp);
      cout << ww[i][j] << " ";
    }
    cout << endl;
  }
  // Now finish with a robustification
  if (robust) {
    // get the diaganol components of the inverted matrix
    vector<double> diag_inv;
    for (i = 0; i < inv_Vmat.size1(); ++i)
      diag_inv.push_back(inv_Vmat(i,i));
    for (unsigned year = 0; year < N_years; ++year) {

    }


  } else {
    vector<vector<double>> temp1,temp2;
    double temp3;
    for (unsigned year = 0; year < N_years; ++year) {
      temp1 = mat_multiply(ww[year], V_invert);
      temp2 = mat_multiply(temp1,ww[year]);
      temp3 = 0.5 / (weights[year] * weights[year]);
      cout << "temp3 = " << temp3 << " temp2 " <<  temp2[0][0] << " current NLL " << score << endl;
      score += (temp3 * temp2[0][0]);
    }
  }



  cout << "--------------------------------------"<< endl;
  cout << "Negative loglikelihood----------------"<< endl;
  cout << "--------------------------------------"<< endl;
  cout << score << endl;



  return score;
}

vector<double> rlogistnorm(vector<double> expprop, double sigma, vector<double> phi, vector<vector<double>> covmat, bool ARMA) {
  vector<double> sim_data;
  sim_data.push_back(0.0);
  return sim_data;
}


// Axuillary fucntion definitions
vector<vector<double> > covmat_logistnorm(double sigma, vector<double> phi,unsigned N_ages , bool sepbysex, bool sexlag, bool ARMA) {
  unsigned n_phi = phi.size();
  vector<vector<double>> covar;
  covar.resize(N_ages, vector<double>(N_ages,0.0));
  if (phi[0] == 0.0) {
    for (unsigned diag = 0; diag < N_ages; ++ diag)
      covar[diag][diag] =sigma * sigma;
  } else {
    vector<double> rho;
    // Get Rho
    cout << "about to enter getRho(), n_phi = " << n_phi << " " << phi[0]<<endl;
    if (n_phi == 1) {
      rho = GetRho(phi[0],N_ages);
    } else if (n_phi == 2 && ARMA) {
      rho = ARMAacf(phi[0],phi[1],N_ages);
    } else {
      rho = ARMAacf(phi,N_ages);
    }

    cout << "Check out rho vector";
    for ( auto num : rho)
      cout << num << " ";
    cout << endl;

    // Simple unsexed example currently
    // Create an identy matrix.
    for (unsigned diag = 0; diag < N_ages; ++ diag) {
      covar[diag][diag] = 1.0 * sigma * sigma;
    }
    for (unsigned diag = 0; diag < N_ages; ++ diag) {
      for (int row = 0; row < N_ages; ++ row) {
        for (int col = 0; col < N_ages; ++ col) {
          if (fabs(row - col) == diag + 1)
            covar[row][col] = rho[diag] * sigma * sigma;
        }
      }
    }
  }
	return covar;
}

// Compute the theoretical autocorrelation function or partial autocorrelation function for an ARMA process.
// The methods used follow Brockwell & Davis (1991, section 3.3). Their equations (3.3.8) are solved for the autocovariances at lags 0, …, max(p, q+1), and the remaining autocorrelations are given by a recursive filter.

vector<double> ARMAacf(double AR, double MA, int nBin) {
  cout << "Beginning method ARMAacf() " << endl;

  // Create and declare all variables that will be used in the function
  unsigned p = 1;
  unsigned q = 1;
  vector<double> AR_coef, MA_coef,final_acf,Cor;
  AR_coef.push_back(AR);
  MA_coef.push_back(MA);

  vector<vector<double> > A, ind;
  vector<double> psi, theta, Acf;
  if (!p && !q)
    cerr << "empty model supplied" << endl;
  unsigned r = fmax(p, q + 1);

  for (unsigned i = 0; i <= (r - p); ++i) {
    AR_coef.push_back(0.0);
    p = r;
  }
  A.resize(p + 1, vector<double>(2 * p + 1, 0.0));
  ind.resize(2 * p + 1, vector<double>(p + 1, 0.0));
  for (unsigned i = 0; i < ind.size(); ++i) {
    for (unsigned j = 0; j < ind[i].size(); ++j) {
      ind[i][0] = i + 1;
      ind[i][1] = (i + 1) + (j + 1) - (i + 1);
    }
  }
  for (unsigned i = 0; i < A.size(); ++i) {
    A[i][i] = 1.0;
    A[i][i + 1] = -AR_coef[0];

  }
  A[0][A.size() - 1] = -AR_coef[1];
  A[A.size() - 1][0] = -AR_coef[1];
  A[A.size() - 1][1] = -AR_coef[0];


  //cout << "size of rhs " << rhs.size() << endl;
  psi.push_back(1.0);
  psi.push_back(AR + MA);
  theta.push_back(1.0);
  theta.push_back(MA);
  for (int i = 0; i <= q; ++i)
    theta.push_back(0.0);

  // Declare Eigen variables
  Eigen::Matrix3d A_eig; // 3f means 3x3 elements of type double
  Eigen::Vector3d rhs;
  // Calculate rhs
  for (unsigned i = 0; i <= q; ++i) {
    double x1 ,x2,tot = 0;
    x1 = psi[0]*theta[i];
    x2 = psi[1]*theta[i + q];
    rhs(i) = Sum({x1,x2});
  }
  rhs(2) = 0.0;

  // Use the eigen library yo solve the inverse of for A with known vector B
  //vector<double> Ind;
  vector<unsigned> seq;
  for (unsigned i = 0; i <= p; ++i) {
    seq.push_back(p - i);
  }


  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p; ++j) {
      cout << ": i = " << i << " j = " << j << " i index = " << seq[i] << " j index = " << seq[j] << " mat value = " << A[seq[i]][seq[j]] << endl;
      A_eig(i,j) = A[seq[i]][seq[j]];
    }
  }

  cout << "Check A mat that we are inverting\n" << A_eig << "\n: rhs = " << rhs << endl;

  // Solve for A_eig given rhs
  Eigen::Vector3d x = A_eig.colPivHouseholderQr().solve(rhs);
  cout << "solution = " << x;

  // Divide the last two elements by the first element.
  cout << "Find the crash" << endl;

  for (unsigned i = 1; i <= 2; ++i) {
    final_acf.push_back(x(i) / x(0));
  }

  cout << "Final Acf" << endl;
  for (auto num : final_acf)
    cout << num << " ";
  cout << endl;

  Cor = RecursiveFilter(AR_coef, nBin, final_acf);

  // Add the initial coeffiecients to the beginning of the series.
  Cor.insert(Cor.begin(), final_acf[1]);
  Cor.insert(Cor.begin(), final_acf[0]);
  // Print results to screen
  vector<double>::const_iterator first = Cor.begin();
  vector<double>::const_iterator last = Cor.begin() + nBin;
  vector<double> newVec(first, last);
  for (auto num : newVec)
    cout << num << " ";
  cout << endl;

  return newVec;

}

vector<double> ARMAacf(vector<double> AR, int nBin) {
  cout << "Beginning method ARMAacf() " << endl;
  // Create and declare all variables that will be used in the function
  unsigned p = AR.size();
  if (p > 2) {
    cerr << "This function has not been coded for more the two AR coeffs." << endl;
  }
  vector<vector<double> > A;
  vector<double> rhs, psi, theta, Acf;
  unsigned r = p;
  A.resize(p + 1, vector<double>(2 * p + 1, 0.0));

  for (unsigned i = 0; i < A.size(); ++i) {
    A[i][i] = 1.0;
    A[i][i + 1] = -AR[0];
    A[i][i + 2] = -AR[1];
  }

  rhs.assign(p + 1, 0.0);
  rhs[0] = 1.0;


  cout << "size of rhs " << rhs.size() << endl;
  //vector<double> Ind;
  vector<unsigned> seq;
  for (unsigned i = 0; i <= p; ++i) {
    seq.push_back(p - i);
  }
  // Create a ublas matrix object for the inversion
  boost::numeric::ublas::matrix<double> A_inv(3,3);
  boost::numeric::ublas::matrix<double> A_ublas(3,3);

  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p ; ++j) {
      if (j == 2)
        A_ublas(i,j) = A[i][j];
      else
        A_ublas(i,j) = A[i][j] + A[i][2 * p  - j];
    }
  }

  for (unsigned i = 0; i <= p; ++i) {
    for (unsigned j = 0; j <= p ; ++j) {
        A_ublas(i,j) =  A_ublas(seq[i],seq[j]);
    }
  }
  // the bodge
  A_ublas(1,2) = 0.0;

  cout << "Check A mat that we are inverti"
      "ng" << endl;
  for(unsigned i = 0; i < A_ublas.size1();++i) {
    for(unsigned k = 0; k < A_ublas.size1();++k) {
      cout << A_ublas(i,k) << " ";
    }
    cout << endl;
  }
  bool inverted = InvertMatrix(A_ublas,A_inv);
  if (!inverted) {
    cerr << "could not invert the matrix" << endl;
  }
  cout << "Check inverted matrix" << endl;
  for(unsigned i = 0; i < A_inv.size1();++i) {
    for(unsigned k = 0; k < A_inv.size1();++k) {
      cout << A_inv(i,k) << " ";
    }
    cout << endl;
  }

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
  xx.assign(nBin - p, 0.0);
  Cor = RecursiveFilter(AR,nBin, final_acf);
  // Add the initial coeffiecients to the beginning of the series.
  Cor.insert(Cor.begin(),final_acf[1]);
  Cor.insert(Cor.begin(),final_acf[0]);
  // Print results to screen
  vector<double>::const_iterator first = Cor.begin();
  vector<double>::const_iterator last = Cor.begin() + nBin;
  vector<double> newVec(first, last);
  for (auto num : newVec)
    cout << num << " ";
  cout << endl;
  return newVec;
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
  for (unsigned i = 1; i < nBins + 1; ++i) {
    if (i == 1) {
      store_vec[i] =   store_vec[i - 1] *ar_coef[0]  + store_vec[i] *  ar_coef[1];
    } else {
      store_vec[i] = store_vec[i - 1] *  ar_coef[0] + store_vec[i - 2] * ar_coef[1];
    }
    cout << "value = " << store_vec[i];
  }
  // remove the first value
  store_vec.erase(store_vec.begin());
  return store_vec;
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
    cout << rho[i - 1] << " ";
  }
  cout <<  "\n\n\n";

  return rho;
}


//////////////////////
// Utility Functions
//////////////////////

// Utilty funciton for summing
double Sum(vector<double> x) {
    cout << "Beginning method ARMAtoMA() " << endl;
  double Total;
    for (unsigned i = 0; i< x.size(); ++i) {
     Total += x[i];
    }
    cout << "leaving method ARMAtoMA() " << endl;
    return Total;
}

double Mean(vector<double> x) {
    cout << "Beginning method ARMAtoMA() " << endl;
  double Total;
    for (unsigned i = 0; i< x.size(); ++i) {
     Total += x[i];
    }
    cout << "leaving method ARMAtoMA() " << endl;
    return Total / (double)x.size();
}

bool all_ones(vector<double> x) {
  for(auto num : x) {
    if (num != 1.0)
      return false;
  }
  return true;

}
vector<double> elem_prod(vector<double> x,vector<double> y) {
  vector<double> result;
  if (x.size() != y.size()) {
    cerr << "cannot apply this function with uneven sized vectors" << endl;
  }
  for (int i = 0; i < x.size(); ++i) {
    result.push_back(x[i] * y[i]);
  }
  return result;
}

// A method for matrix multiplication
vector<vector<double> > mat_multiply(vector<vector<double> >& x,vector<vector<double> >& y) {
  cout << "Entering mat_multiply() method" <<endl;
  vector<vector<double>> result;
  result.resize(x.size(), vector<double>(x.size()));
  cout << "rows of x = " <<  x.size() << " rows of y = " <<  y.size() << endl;
   for (int i = 0; i < x.size(); i++) {
    for (int k = 0; k < x[i].size(); k++) {
     for (int j = 0; j < y.size(); j++) { // swapped order
       result[i][j] += x[i][k] * y[k][j];
     }
    }
   }
   return result;
}

// A method for matrix multiplication
vector<vector<double> > mat_multiply(vector<double> & x,vector<vector<double> >& y) {
  vector<vector<double>> result;
  result.resize(1, vector<double>(x.size()));
   for (int i = 0; i < y.size(); i++) {
    for (int k = 0; k < y[i].size(); k++) {
       result[0][k] += x[i] * y[i][k];
    }
   }
   return result;
}
// A method for matrix multiplication
vector<vector<double> > mat_multiply(vector<vector<double> >& x,vector<double> & y) {
  vector<vector<double>> result;
  result.resize(x.size(), vector<double>(1));
   for (int i = 0; i < x.size(); i++) {
    for (int k = 0; k < x[i].size(); k++) {
       result[i][0] += y[k] * x[i][k];
    }
   }
   return result;
}
vector<vector<double> > t(vector<vector<double> >& x) {
  vector<vector<double>> result;
  result.resize(x[0].size(), vector<double>(x.size()));
  for (int i = 0; i < x.size(); i++) {
    for (int j = 0; j < x[i].size(); j++) {
      result[j][i] = x[i][j];
    }
  }
  return result;
}

vector<vector<double>> log_mat(vector<vector<double> >& x) {
  vector<vector<double>> result;
  result.resize(x.size(), vector<double>(x[0].size()));
  for (unsigned i = 0; i < x.size(); ++i) {
    for (unsigned j = 0; j < x[i].size(); ++j)
      result[i][j] = log(x[i][j]);
  }
  return result;
}

double Sum_mat(vector<vector<double> >& x) {
  double Tot = 0.0;
  for (unsigned i = 0; i < x.size(); ++i) {
    for (unsigned j = 0; j < x[i].size(); ++j)
      Tot += x[i][j];
  }
  return Tot;
}

// solve for the inverse of n x n matrix
bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse) {
  typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  boost::numeric::ublas::matrix<double> A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = boost::numeric::ublas::lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(boost::numeric::ublas::identity_matrix<double> (A.size1()));

  // backsubstitute to get the inverse
  boost::numeric::ublas::lu_substitute(A, pm, inverse);

  return true;
}

// Calculate the determinant
double det_fast(const boost::numeric::ublas::matrix<double>& matrix) {
    // create a working copy of the input
  boost::numeric::ublas::matrix<double> mLu(matrix);
  boost::numeric::ublas::permutation_matrix<std::size_t> pivots(matrix.size1());

    auto isSingular = boost::numeric::ublas::lu_factorize(mLu, pivots);
    if (isSingular)
        return static_cast<double>(0);

    double det = static_cast<double>(1);
    for (size_t i = 0; i < pivots.size(); ++i)
    {
        if (pivots(i) != i)
            det *= static_cast<double>(-1);

        det *= mLu(i, i);
    }

    return det;
}

