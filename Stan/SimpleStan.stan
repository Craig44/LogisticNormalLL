/*
 * @file AgeStructuredModel_1.0.stan
 * @author C. Marsh
 *
 * A stan module for a simple age structured dynamics model
 */
 
 
functions {
  // logistic ogive function
  real gmean(row_vector x) {
    real tot = 1.0;
    for (i in 1:num_elements(x))
      tot *= x[i];
    return tot^(1.0/num_elements(x));
  }
  /* compute the cholesky factor of an AR1 correlation matrix
   * @author paul.buerkner@gmail.com
   * Args: 
   *   ar: AR1 autocorrelation 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows matrix 
   */ 
   matrix cholesky_cor_ar1(real ar, int nrows) { 
     matrix[nrows, nrows] mat; 
     vector[nrows - 1] gamma; 
     mat = diag_matrix(rep_vector(1, nrows)); 
     for (i in 2:nrows) { 
       gamma[i - 1] = pow(ar, i - 1); 
       for (j in 1:(i - 1)) { 
         mat[i, j] = gamma[i - j]; 
         mat[j, i] = gamma[i - j]; 
       } 
     } 
     return cholesky_decompose(mat); 
   }
  // Logisitic likelihood negative loglikelihood
  real calculate_logistic_lpdf(matrix Observations, matrix expectations, matrix covariance, vector sample_size_by_year) {
    int A = cols(expectations);
    int Y = rows(expectations);
    matrix[A - 1, A - 1] v_mat;
    matrix[A - 1, A - 1] inv_v_mat;
    matrix[A - 1, A] k_mat;
    matrix[A, A - 1] t_k_mat;
  
    matrix[Y, A-1]  ww;
    vector[Y] weights_by_year;  // Sample sizes to 
    real log_det;
    real temp_val = 0.0;
    real log_total_obs = 0.0;
    
    real neg_ll;
    // initialisa
    k_mat = rep_matrix(0.0, A - 1, A);
    v_mat = rep_matrix(0.0, A - 1, A - 1);
    ww = rep_matrix(0.0, Y, A - 1);
    k_mat[,A] = rep_vector(-1.0, A - 1);
    for (i in 1:(A-1)) 
      k_mat[i,i] = 1.0;
    t_k_mat = covariance * (k_mat');
    v_mat = k_mat * t_k_mat;
    inv_v_mat = inverse(v_mat);
  
    log_det = log_determinant(v_mat);
    log_total_obs = sum(log(Observations));
    neg_ll = 0.5 * Y * (A - 1)* log(2.0*pi()) + log_total_obs + 0.5 * Y * log_det;
    weights_by_year = sqrt(mean(sample_size_by_year) ./ sample_size_by_year);
    neg_ll += (A - 1) * sum(log(weights_by_year));
    for (y in 1:Y) 
      ww[y,] = log(Observations[y, 1:(A-1)] ./ Observations[y, A]) - log(expectations[y, 1:(A-1)] ./ expectations[y, A]);
    
    for (y in 1:Y) 
      neg_ll += sum((0.5/(weights_by_year[y]*weights_by_year[y])) * (ww[y,] * inv_v_mat) .* ww[y,]);
    return neg_ll;
  }
  
  matrix standardised_residuals(matrix Observations, matrix expectations, matrix covariance, vector sample_size_by_year, int centered) {
    int A = cols(expectations);
    int Y = rows(expectations);
    matrix[A - 1, A - 1] v_mat;
    matrix[A - 1, A - 1] inv_v_mat;
    matrix[A - 1, A] k_mat;
    matrix[A - 1, A] FF;
    matrix[A, A - 1] t_k_mat;
    matrix[A-1, A-1] HH;
    matrix[Y, A-1]  ww;
    matrix[A, A-1] FHinv;
    matrix[A, A] Gam;
    matrix[Y, A] sdmat;
    matrix[Y, A - 1] sdmat_non_centered;
    matrix[Y, A] sres;
    matrix[Y, A - 1] sres_non_centered;

    vector[Y] gmean_obs;  // Sample sizes to 
    vector[Y] gmean_exp;
    vector[Y] weights_by_year;  // Sample sizes to 
    real log_det;
    real temp_val = 0.0;
    real log_total_obs = 0.0;
    
    // initialisa
    k_mat = rep_matrix(0.0, A - 1, A);
    v_mat = rep_matrix(0.0, A - 1, A - 1);
    ww = rep_matrix(0.0, Y, A - 1);
    weights_by_year = sqrt(mean(sample_size_by_year) ./ sample_size_by_year);

    k_mat[,A] = rep_vector(-1.0, A - 1);
    for (i in 1:(A-1)) 
      k_mat[i,i] = 1.0;

    
    t_k_mat = covariance * (k_mat');
    v_mat = k_mat * t_k_mat;
    inv_v_mat = inverse(v_mat);
    
    if(centered) {
        HH = rep_matrix(1.0, A - 1, A - 1);
        for (i in 1:cols(HH))
          HH[i,i] = 2.0;
          
        FF  = k_mat;
        FF[,A] = rep_vector(1.0, A - 1);
        FHinv = FF' * inverse(HH);
        Gam = FHinv * (v_mat * FHinv');
        for(y in 1:Y) 
          sdmat[y,] = (sqrt(diagonal(Gam)) * weights_by_year[y])';
        for(y in 1:Y)  {
          gmean_obs[y] = gmean(Observations[y,]);
          gmean_exp[y] = gmean(expectations[y,]);
          sres[y, ] = (log(Observations[y,] / gmean_obs[y]) - log(expectations[y,] / gmean_exp[y])) ./ sdmat[y,];
        }
        return sres;
    } else {
      for(y in 1:Y) 
        sdmat_non_centered[y,] = (sqrt(diagonal(v_mat)) * weights_by_year[y])';
      for(y in 1:Y)  {
        sres_non_centered[y, ] = (log(Observations[y,1:(A - 1)] / Observations[y,A]) - log(expectations[y,1:(A - 1)] / expectations[y,A])) ./ sdmat_non_centered[y,];
      }
      return sres_non_centered;
    }
    
    // Should have exited before here, but compiler will complain.
    return sres;
  }

  matrix ar1_covar(real sigma, real rho, int bins) {
    matrix[bins,bins] covariance_matrix;
    matrix[bins,bins] cholesky_cor = cholesky_cor_ar1(rho, bins);
    covariance_matrix = diag_pre_multiply(rep_vector(sigma,bins), cholesky_cor) * diag_pre_multiply(rep_vector(sigma,bins), cholesky_cor)';
    return covariance_matrix;
  }

}
data {
  // Model dimensions
  int<lower=1> Y; // number of years
  int<lower=1> A; // number of ages
  // Observational inputs
  matrix[Y,A] Observations; // Observations
  matrix[Y,A] expectations; // Observations
  matrix[A,A] covariance;
  vector[Y] sample_size_by_year;  // Sample sizes to 
  int<lower=0, upper = 1> centered; // if = 1 - standardised residuals are based on the MVN variate Z, where Zi = log(Xi/GM(X)); otherwise (=0) they are based on the MVN variate Y, where Yi = log(Xi/XD), where D = ncol(obsprop)
  real rho;
  real sigma;

}

// This section is equivalent to dobuild and do validate, for operations that don't need to be
// re executed they should go here. For example, truncating observations etc
transformed data {
  
  
}

/*
 * The estimated parameters
*/
parameters {
  real<lower = 0, upper = 1> p;
}

/*
 * Transformed Data
 * The purpose of this section of code is for ...
*/

transformed parameters{
  // Define variables
  real neg_ll = 0.0;
  // Build covariance matrix
  matrix[A,A] covariance_matrix;
  covariance_matrix = ar1_covar(sigma, rho, A);
  // Calculates standardised residuals for composition data
  neg_ll = calculate_logistic_lpdf(Observations | expectations, covariance_matrix, sample_size_by_year);
}

model {
  // turn it into a positive because stan does maximum likelihood not neg ll
  target += (1.0);
}

generated quantities {
  matrix[Y,A] sim_bservations;
  for (y in 1:Y) {
   sim_bservations[y,] = softmax(multi_normal_rng(expectations[y,],covariance_matrix))';
  }
}
