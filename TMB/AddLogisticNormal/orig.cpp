#include <TMB.hpp>
#include "Logistic_normal_likelihood.hpp"

// Square a quantity
template <class Type> 
Type square(Type x){
  return x*x;
}

// Logistic transform
template <class Type> 
Type logistic_transform(Type x) {
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}


// Zerofun functions
template <class Type> 
Type ZeroFun(Type x, Type delta) {
  if (x >= delta)
    return x;
  
  return delta / (2.0 - (x / delta));
}


// logistic ogive function
template <class Type> 
vector<Type> logistic_ogive(vector<Type> ages, Type sel_50, Type sel_95) {
  std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + pow(Type(19.0), (sel_50 - ages[age]) / sel_95));
  }
  return logis;
}

// Weight at age
template <class Type> 
vector<Type> mean_weight_at_age(vector<Type> L_a, Type a, Type b) {
  std::cout << "mean_weight\n";
  int n_ages = L_a.size();
  vector<Type> W_a(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    W_a(age) = a * pow(L_a(age), b);
  }
  return W_a;
}
 
 // Length at age
 template <class Type> 
 vector<Type> VonBertalanffy(vector<Type>& ages, Type& L_inf, Type& k, Type& t0) {
   std::cout << "VB\n";
   int n_ages = ages.size();
   vector<Type> L_a(n_ages);
   for (int age = 0;  age < n_ages; ++age) {
     L_a(age) = L_inf * (Type(1) - exp(-k*(ages(age) -t0)));
   }
   return L_a;
 }   

// Beverton-Holt SR relationship function
template<class Type>
Type BevertonHolt(Type SSB, Type B0, Type h) {
  Type ssb_ratio = SSB / B0;
  Type part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
  return (ssb_ratio / part_2);
}

/*Likelihood for CASAL comparisons*/

// objective calculations, from CASAL
// returns negative log-likelihood.
template <class Type> 
Type obj_multinomal_casal_log(array<Type> obs, matrix<Type> Exp, vector<Type> error) {
  std::cout << "obj_multinomal_casal_log\n";
  int Y = Exp.rows();
  int A = Exp.cols();
  Type obj = 0.0;
  for (int y = 0; y < Y; ++y) {
    obj = obj - lgamma(error[y] + 1.0); // lgamma(n + 1) = lfactorial(n)
    //obj += lgamma(error[y] + 1.0); // lgamma(n + 1) = lfactorial(n)
    for (int age = 0;  age < A; ++age) {
      obj = obj + lgamma(error(y) * obs(y,age) + 1.0) - error(y) * obs(y,age) * log(Exp(y,age));
    }
  }
  return obj;
}

// return pearsons residuals
template <class Type> 
array<Type> pearsons_resids_multinomial(array<Type> obs, matrix<Type> Exp, vector<Type> error) {
  int Y = Exp.rows();
  int A = Exp.cols();
  array<Type> resids(obs);
  for (int y = 0; y < Y; ++y) {
    for (int age = 0;  age < A; ++age) {
      resids(y,age) = (obs(y,age) - Exp(y,age)) / sqrt((ZeroFun(Exp(y,age), Type(0.000001)) * (1 - ZeroFun(Exp(y,age), Type(0.000001)))) / error(y));
    }
  }
  return resids;
} 

// negative log likelihood for the lognormal distribution
template <class Type> 
Type obj_lognomal_like(vector<Type> observed, vector<Type> expected, vector<Type> error) {
  std::cout << "obj_lognomal_like\n";
  int n_years = observed.rows();
  vector<Type> sigma(n_years);
  vector<Type> score(n_years);
  Type neg_ll = 0.0;
  for (int y = 0; y < n_years; ++y) {
    sigma(y) = sqrt(log(1.0 + error(y) * error(y)));
    score(y) = log(observed(y) / expected(y))/ sigma(y) + 0.5 * sigma(y);
    neg_ll += log(sigma(y)) + 0.5 * score(y) * score(y);
  }
  return neg_ll;
}  
// return normalised residuals
template <class Type> 
vector<Type> norm_resids_lognormal(vector<Type> observed, vector<Type> expected, vector<Type> error) {
  int n_years = observed.size();
  vector<Type> sigma(n_years);
  vector<Type> resids(n_years);
  for (int y = 0; y < n_years; ++y) {
    sigma(y) = sqrt(log(1.0 + error(y) * error(y)));
    resids(y) = (log(observed(y) / expected(y)) + 0.5 * sigma(y) * sigma(y)) / sigma(y);
  }
  return resids;
} 
// Casal's lognormal prior calculation
template <class Type> 
Type lnorm_prior(Type X, Type Expectation, Type sigma) {
  Type result = log(X) + log(sigma) + 0.5 * pow(log(X/ Expectation) / sigma + sigma * 0.5,2);
  return result;
}




/*
* Simplex functionality
*/
// simplex_jacobian - calculate the jacobian for estiamting the logistic simplex rather than the simplex
template <class Type> 
Type simplex_jacobian(vector<Type> zk, vector<Type> X_orig) {
  std::cout << "Simplex Jacobian\n";
  Type Jacobian = 1.0;
  for (int k = 0; k < (zk.size() - 1); ++k) { // because I store zk in a matrix with xk in simplex_restore() the zk vector has one more null/empty element than it should have so we don't iterate over this.
    if (k == 0) {
      Jacobian *= (zk[k] * (1 - zk[k]));  
    } else {
      Jacobian *= (zk[k] * (1 - zk[k])) * (1 - sum(vector<Type>(X_orig.segment(0,k))));
    }
  }
  return Jacobian;
}

// This class takes the transformed (logistic simplex) values, and their original unit vector and gives the
// recent update value.
template <class Type>
matrix<Type> Simplex_restore(vector<Type> logit_unit_vec, vector<Type> last_unit_vec) {
  std::cout << "Simplex_restore\n";
  int n_pars = last_unit_vec.size();
  matrix<Type> matrix_return(2,n_pars);
  for (int k = 0; k < (logit_unit_vec.size()); ++k) {
    matrix_return(0,k) = invlogit(logit_unit_vec[k] + log(Type(1) / ((Type(logit_unit_vec.size()) + Type(1)) - Type(k + 1))));
    
    if (k == 0)
      matrix_return(1,k) = matrix_return(0,k);
    else 
      matrix_return(1,k) = (Type(1) - sum(vector<Type>(last_unit_vec.segment(0,k)))) * matrix_return(0,k);
    //std::cerr << "k = " << k + 1 << " zk = " << matrix_return(0,k) << " xk = " << matrix_return(1,k) << std::endl; // Output warning
  }
  //std::cerr << vector<Type>(matrix_return.row(1)) << std::endl;
  matrix_return(1,last_unit_vec.size() - 1) = 1 - sum(vector<Type>(matrix_return.row(1)));
  //std::cerr << vector<Type>(matrix_return.row(1)) << std::endl;
  
  return matrix_return;
}


template<class Type>
Type objective_function<Type>::operator() () {
  // Model Dimensions
  DATA_VECTOR(years);
  DATA_VECTOR(ages); // number of ages
  int Y = years.size();
  int A = ages.size();
  // Observations  
  DATA_ARRAY(fishery_at_age_obs);
  DATA_VECTOR(fishery_at_age_error);
  DATA_ARRAY(survey_at_age_obs);
  DATA_VECTOR(survey_at_age_error);
  DATA_VECTOR(survey_biomass_obs);
  DATA_VECTOR(survey_biomass_error);
  DATA_VECTOR(fishery_years);
  DATA_VECTOR(survey_years);
  DATA_MATRIX(ageing_error);
  DATA_MATRIX(fishery_at_age_exp_LN_test);
  DATA_VECTOR(fishery_at_age_error_LN_test);
  
  // Biological parameters.
  DATA_SCALAR(h);// Steepnees
  DATA_SCALAR(a);// a in the mean weight calculation
  DATA_SCALAR(b);// b in the mean weight calculation
  DATA_SCALAR(L_inf);// L_inf in the Von Bertallanffy equation
  DATA_SCALAR(k);// k in the Von Bertallanffy equation
  DATA_SCALAR(t0);// t0 in the Von Bertallanffy equation
  DATA_SCALAR(mat_a50);// Maturity Parameters
  DATA_SCALAR(mat_ato95);// Maturity Parameters
  DATA_SCALAR(M);
  
  // Other flags and variables
  DATA_INTEGER(ycs_prior_applies_to_constrained);// 0 = no,1 = yes priors apply to constrained or unconstrained parameters.
  DATA_INTEGER(standardise_ycs); // // 0 = no,1 = yes apply standardisation on YCS
  DATA_INTEGER(simplex);  // // 0 = no,1 = yes apply the simplex transformation on YCS
  DATA_INTEGER(deviations); // // 0 = no,1 = yes YCS or devs, if devs it is assumed that they YCS vector are constrained sum = 0
  DATA_INTEGER(use_logistic_normal); //should comp data be modelled by the logistic normal
  DATA_INTEGER(ARMA); // Specifies a certain correlation structure in the covariance.
  
  DATA_VECTOR(untransformed_values); // The starting values if doing the simplex transformation
  
  // Exploitation crap
  DATA_SCALAR(u_max);
  DATA_SCALAR(Catch_penalty_multiplier);//
  DATA_SCALAR(proportion_mortality_spawning);
  DATA_SCALAR(proportion_mortality_survey);
  DATA_VECTOR(catches);
  
  // Priors/penalties
  DATA_INTEGER(catch_penalty_log_space); // penalty flag
  DATA_INTEGER(apply_priors); // Do we want to apply priors/penalties to YCS and or Q's

  DATA_SCALAR(mu_q);
  DATA_SCALAR(sigma_q);
  
  // One step predictions
  DATA_ARRAY_INDICATOR(survey_age_pred, survey_at_age_obs);
  DATA_ARRAY_INDICATOR(fishery_age_pred, fishery_at_age_obs);
  DATA_VECTOR_INDICATOR(survey_biomass_pred, survey_biomass_obs);  // For one-step predictions
  
  
  // Estimated parameters
  PARAMETER(log_R0);
  PARAMETER(q);
  PARAMETER(s_a50);
  PARAMETER(s_ato95);
  PARAMETER(f_a50);
  PARAMETER(f_ato95);
  PARAMETER(log_sigma_r);
  PARAMETER_VECTOR(YCS);
    // Logistic normal values
  PARAMETER(log_norm_sigma);
  PARAMETER_VECTOR(log_norm_phi);
  /*
   * Transformations
   */
  vector<Type> true_ycs;
  Type R0 = exp(log_R0);
  
  Type sigma_r = exp(log_sigma_r);
  
  if ((standardise_ycs == 1) && (simplex == 1) ) {
    error("You cannot supply both simplex and standardise ycs");
  }
  
  if ((deviations == 1) && (simplex == 1) ) {
    error("You cannot supply both simplex and deviations");
  }  
  
  true_ycs = YCS;
  
  matrix<Type> untransformed_simplex;  
  Type jacobian_from_simplex = 1.0;
  if (simplex == 1) {
    if ((YCS.size() + 1) != Y) 
      error("you must supply one less simplex than year, as it is constrained");
    
    untransformed_simplex = Simplex_restore(YCS, untransformed_values);
    //untransformed_values = vector<Type>(untransformed_simplex.row(1));
    jacobian_from_simplex = simplex_jacobian(vector<Type>(untransformed_simplex.row(0)), untransformed_values);
    untransformed_values = vector<Type>(untransformed_simplex.row(1));
    
    true_ycs = untransformed_values * Type(Y);
    std::cerr << true_ycs << std::endl;
    std::cerr << true_ycs.size() << std::endl;
    //return 0;
  }
  vector<Type> transformed_deviations(Y);
  if (deviations == 1) {
    std::cerr << "ycs = " << (YCS.size() + 1) << " Y = " << Y << std::endl;
    if ((YCS.size() + 1) != Y) 
      error("you must supply one less deviation than year, as it is constrained");
    for (int i = 0; i < (YCS.size()); ++i) {
      transformed_deviations(i) = YCS(i);
    }
    transformed_deviations(YCS.size()) = -sum(YCS);
    true_ycs = exp(transformed_deviations - 0.5 * sigma_r * sigma_r);

    std::cerr << "29th dev = " << transformed_deviations(YCS.size() - 1) << std::endl;
    
    
    std::cerr << "sum of devs = " << sum(transformed_deviations) << std::endl;
    std::cerr << "size of devs = " << transformed_deviations.size() << std::endl;
    std::cerr << "devs = " << transformed_deviations << std::endl;
    std::cerr << "size of y = " << true_ycs.size() << std::endl;
    std::cerr << "true_ycs = " << true_ycs << std::endl;
    
    
    if (fabs(sum(transformed_deviations)) > 0.0001)
      error("deviations are constrained to sum to 0 but this isn't the case please address");
  }
  
  if (standardise_ycs == 1) {
    true_ycs = YCS / YCS.mean();
  } 
  
  
  /*
   * Big list of sanity checks that we would want to add to check before going anyfurtur
   */
  if (ageing_error.rows() != ageing_error.cols())
    error("ageing error needs to be symetric, ie rows = cols");
  if (ageing_error.rows() != A)
    error("ageing error needs to have rows and cols = number of ages this is not the case");
    
  
  
  
  
  
  
  // Define containers that will be needed over the model life.
  // Matrix
  matrix<Type> numbers_at_age(Y+1,A);
  matrix<Type> pre_survey_age_expectations(Y, A);
  matrix<Type> survey_age_expectations(survey_years.size(), A);
  matrix<Type> fishery_age_expectations(fishery_years.size(), A);
  

  numbers_at_age.setZero();
  pre_survey_age_expectations.setZero();
  fishery_age_expectations.setZero();
  survey_age_expectations.setZero();
  // Vectors used
  vector<Type> temp_partition(A); // vector of ages
  vector<Type> recruits(Y);
  vector<Type> SSB(Y);
  vector<Type> maturity_at_age(A); // maturity schedule
  vector<Type> fish_select_at_age(A); // Fishery selectivity
  vector<Type> survey_select_at_age(A); // survey selectivity
  vector<Type> length_at_age(A); // length at age
  vector<Type> weight_at_age(A); // weight at age
  vector<Type> pre_SSB(Y);
  vector<Type> pre_survey_biomass(Y); // I don't do a check on pre survey calculations wheter an observation exists thus this needs to exist over all years
  vector<Type> vulnerable(Y);
  vector<Type> exploitation(Y);
  vector<Type> actual_catches(Y);
  vector<Type> survey_biomass_expectations(survey_years.size());

  // SCalars used
  Type B0;
  Type pre_B0;
  Type u_obs;
  Type neg_ll_fishery_age = 0;
  Type neg_ll_survey_age = 0;
  Type neg_ll_survey_bio = 0;
  Type penalty = 0;
  

  // Begin calcualtions

  std::cout << "R0 q s_a50 s_ato95 f_a50 f_ato95 YCS\n";
  std::cout << R0 << " " << q << " " << s_a50 << " " << s_ato95 << " " << f_a50 << " " << f_ato95 << " " << YCS << "\n";

  // Do preliminary calculations
  length_at_age = VonBertalanffy(ages, L_inf, k, t0);
  weight_at_age = mean_weight_at_age(length_at_age, a, b);
  fish_select_at_age = logistic_ogive(ages, f_a50, f_ato95);
  maturity_at_age = logistic_ogive(ages, mat_a50, mat_ato95);
  survey_select_at_age = logistic_ogive(ages, s_a50, s_ato95);

  /* 
   *  ======================
   *  Initialise Partition
   *  ======================
   */   
  //
  numbers_at_age(0,0) = R0;
  for(int age = 1; age < (A - 1); ++age) {
    numbers_at_age(0, age) = numbers_at_age(0, age - 1) * exp(-M);
  }
  numbers_at_age(0,A - 1) =  numbers_at_age(0,A - 2) * exp(-M) / (1 - exp(-M));
  pre_B0 = sum(vector<Type>(vector<Type>(numbers_at_age.row(0)) * maturity_at_age * weight_at_age));
  temp_partition = numbers_at_age.row(0);
  numbers_at_age.row(0) = temp_partition * exp(-M);
  vector<Type> initial = numbers_at_age.row(0);
  B0 = pre_B0 + (sum(vector<Type>(vector<Type>(numbers_at_age.row(0)) * maturity_at_age * weight_at_age)) - pre_B0) * proportion_mortality_spawning;
  
  /* 
   *  ======================
   *  run model for n years
   *  ======================
   */
  for (int t = 0; t < Y; ++t) {  
    /////////////////
    // Ageing
    /////////////////
    temp_partition = numbers_at_age.row(t);
    for(int age = 0; age < (A - 2); ++age) 
      numbers_at_age(t + 1,age + 1) = numbers_at_age(t, age);
    // plus group
    numbers_at_age(t + 1, A - 1) = temp_partition(A - 2) + temp_partition(A - 1);
    /////////////////
    // Recruitment
    /////////////////    
    if (t == 0) {
      recruits(t) = R0 * true_ycs(t);
    } else {
      recruits(t) = R0 * BevertonHolt(SSB(t-1), B0, h) * true_ycs(t);
    }
    numbers_at_age(t + 1, 0) = recruits(t);
    
    //////////////////
    //  Pre-Execute
    //////////////////
    pre_SSB(t) = sum(vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) *  maturity_at_age *  weight_at_age));
    pre_survey_biomass(t) = sum(vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) *  survey_select_at_age *  weight_at_age));
    pre_survey_age_expectations.row(t) = vector<Type>(numbers_at_age.row(t + 1)) *   survey_select_at_age;
    
    ///////////////////
    // Mortality
    ///////////////////
    vulnerable[t] = sum(vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) * exp(-0.5 * M) * fish_select_at_age * weight_at_age));
    // Calculate u_obs for each fishery, this is defined as the maximum proportion of fish taken from any element of the partition
    u_obs = max(vector<Type>((catches(t) / vulnerable(t)) * fish_select_at_age));
    
    // Uobs is just for reporting and comparing with u_max
    if (u_obs > u_max) {
      exploitation(t) = (catches(t) / vulnerable(t)) * (u_max / u_obs);
      //flag_catch_penalty = 1.0;
      actual_catches(t) = exploitation(t) * vulnerable(t);
      u_obs = u_max;
      if (catch_penalty_log_space == 1) {
        penalty += pow(log(catches(t)) - log(actual_catches(t)),2); 
      } else {
        penalty += pow(catches(t) - actual_catches(t),2); 
      }
    } else {
      exploitation(t) = u_obs;
      actual_catches(t) = catches(t);
      u_obs = catches(t) / vulnerable(t);
    }
    
    ///////////////////
    // Calcualte Fishery Observation halfway through the Mortality process
    ///////////////////
    for (int i = 0; i < fishery_years.size();++i) {
      if (fishery_years(i) == years(t)) {
        fishery_age_expectations.row(i) = vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) * u_obs * fish_select_at_age * exp(-M * 0.5));
        // TODO add ageing error here
        vector<Type> numbers_at_age_with_error(A);
        temp_partition = fishery_age_expectations.row(i);
        
        for (int j = 0; j < A; ++j) {
          for (int p = 0; p < A; ++p) {
            numbers_at_age_with_error[p] += temp_partition(j) * ageing_error(j,p);
          }
        }
        temp_partition = numbers_at_age_with_error;
        fishery_age_expectations.row(i) = temp_partition / sum(temp_partition);
      }
    }
    
    // apply the mortality
    temp_partition = numbers_at_age.row(t + 1);
    numbers_at_age.row(t + 1) = vector<Type>(temp_partition * exp(-M) * (1 - u_obs * fish_select_at_age));
    
    
    // Calculate SSB which is an approximation through the annual time step
    SSB(t) = pre_SSB(t) + (sum(vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) * maturity_at_age * weight_at_age)) - pre_SSB(t)) * proportion_mortality_spawning;
 
    for (int i = 0; i < survey_years.size();++i) {
      if (survey_years(i) == years(t)) {
        // Calculate Survey expectations
        survey_biomass_expectations(i) = pre_survey_biomass(t) + (sum(vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) * survey_select_at_age * weight_at_age)) - pre_survey_biomass(t)) * proportion_mortality_survey;
        survey_biomass_expectations(i) *= q;
        // This is actually numbers at age until we get to the next section
        survey_age_expectations.row(i) = vector<Type>(pre_survey_age_expectations.row(t)) + ((vector<Type>(vector<Type>(numbers_at_age.row(t + 1)) * survey_select_at_age) - vector<Type>(pre_survey_age_expectations.row(t))) * proportion_mortality_survey);
        // TODO add ageing error here
        vector<Type> numbers_at_age_with_error(A);
        temp_partition = survey_age_expectations.row(i);
        
        for (int j = 0; j < A; ++j) {
          for (int p = 0; p < A; ++p) {
            numbers_at_age_with_error[p] += temp_partition(j) * ageing_error(j,p);
          }
        }
        temp_partition = numbers_at_age_with_error;
        survey_age_expectations.row(i) = temp_partition / sum(temp_partition);
      }
    }

  }
  /*
   * Additive Logistic Normal Likelihood
   */
  if (use_logistic_normal == 1) {
    /*
    int N_bins = fishery_at_age_exp_LN_test.cols();
    matrix<Type> covar;
    covar = covariance_logistic(exp(log_norm_sigma),exp(log_norm_phi),N_bins,ARMA);
    matrix<Type> k_mat(N_bins - 1,N_bins);
    k_mat.setZero();
    matrix<Type> V_mat(N_bins - 1,N_bins);
    V_mat.setZero();
    matrix<Type> ww(N_bins,N_bins - 1);
    ww.setZero();
    for (int i = 0; i < (N_bins - 1); ++i) {
      k_mat(i,i) = 1.0;
      k_mat(i, A - 1) = -1.0;
    }
    matrix<Type> t_kmat = covar * k_mat.transpose();
    
    V_mat = k_mat * t_kmat;
    
    // Calculate log determinant and inverse
    Type log_det_please;
    Type log_det = atomic::logdet(V_mat);
    REPORT(log_det);
    
    matrix<Type> V_mat_inv = atomic::matinvpd(V_mat,log_det_please);
    REPORT(V_mat_inv);
    REPORT(log_det_please);
    
    Type log_obs_tot = fishery_at_age_obs.log().sum();
    Type neg_ll_LN = 0.5 * Y * Type(A - 1) * log(2*3.1415926535) + log_obs_tot + 0.5 * Y * log_det_please;
    REPORT(neg_ll_LN);
    REPORT(log_obs_tot);
    REPORT(k_mat);
    REPORT(V_mat);
    matrix<Type> temp = k_mat.transpose();
    REPORT(temp);
    
    */
    neg_ll_fishery_age = NLLlogistnorm(fishery_at_age_obs, fishery_at_age_exp_LN_test, fishery_at_age_error_LN_test,exp(log_norm_sigma),exp(log_norm_phi),ages,ARMA);
    //REPORT(covar);
    //neg_ll_survey_age = NLLlogistnorm(survey_at_age_obs, survey_age_expectations, survey_at_age_error,exp(log_norm_sigma),exp(log_norm_phi),ages,ARMA);
    // ----------------- Test code for covariance creation -------------------------//
  } else {
    neg_ll_survey_age = obj_multinomal_casal_log(survey_at_age_obs, survey_age_expectations, survey_at_age_error);
    neg_ll_fishery_age = obj_multinomal_casal_log(fishery_at_age_obs, fishery_age_expectations, fishery_at_age_error);
  }
  neg_ll_survey_bio = obj_lognomal_like(survey_biomass_obs, survey_biomass_expectations, survey_biomass_error);
  
  // Calculate residuals
  fishery_age_pred = pearsons_resids_multinomial(fishery_at_age_obs, fishery_age_expectations, fishery_at_age_error);
  survey_age_pred = pearsons_resids_multinomial(survey_at_age_obs, survey_age_expectations, survey_at_age_error);
  survey_biomass_pred = norm_resids_lognormal(survey_biomass_obs, survey_biomass_expectations, survey_biomass_error);
  /* 
   *  ===============
   *  Report section
   *  ===============
   */   
  REPORT(ages);
  REPORT(recruits);
  REPORT(vulnerable);
  REPORT(exploitation);
  REPORT(SSB);
  ADREPORT(SSB)
  REPORT(numbers_at_age);
  REPORT(B0);
  REPORT(survey_select_at_age);
  REPORT(fish_select_at_age);
  REPORT(fish_select_at_age);
  REPORT(weight_at_age);
  REPORT(true_ycs);
  REPORT(survey_biomass_expectations);
  REPORT(survey_age_expectations);
  REPORT(fishery_age_expectations);
  
  REPORT(neg_ll_fishery_age);
  REPORT(neg_ll_survey_bio);
  REPORT(neg_ll_survey_age);
  REPORT(jacobian_from_simplex);
  REPORT(untransformed_values);
  /* 
   *  ===============
   *  Simulate section
   *  ===============
   */   

  Type neg_ll = neg_ll_fishery_age + neg_ll_survey_bio + neg_ll_survey_age + (penalty * Catch_penalty_multiplier) - log(jacobian_from_simplex);
  REPORT(neg_ll);
  if (apply_priors) {
    Type q_prior = lnorm_prior(q, Type(mu_q), Type(sigma_q));
    REPORT(q_prior)
    neg_ll += q_prior;
    Type ycs_prior = 0;
    if (ycs_prior_applies_to_constrained == 1) {
      for (int i = 0; i < true_ycs.size(); ++i) {  
        ycs_prior += lnorm_prior(true_ycs(i), Type(1.0), Type(sigma_r));
      }
    } else if (deviations == 1){
      for (int i = 0; i < transformed_deviations.size(); ++i) {  
        ycs_prior = -dnorm(transformed_deviations(i), Type(0.0), Type(sigma_r), true);
      }
    } else {
      for (int i = 0; i < YCS.size(); ++i) {  
        ycs_prior += lnorm_prior(YCS(i), Type(1.0), Type(sigma_r));
      }
    }

    REPORT(ycs_prior)
    neg_ll += ycs_prior;
    
  }
  return neg_ll;
}

//gdbsource("my_project.R")
