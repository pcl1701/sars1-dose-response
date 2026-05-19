functions {
  
  /**
   * Compute viral load
   *
   * @param t time to compute viral load
   * @param g exponential growth rate
   * @param d exponential decay rate
   * @param tp time of peak viral load
   * @param logN0 initial log viral load
   * @return expected log viral load value at time t
   */
  real log_viral_load(real t, real g, real d, real tp, real logN0) {
    
    real result = 0;
    
    // growth phase
    if (t < tp) {
      result = logN0 + g * t;
    }
    // decay phase
    else {
      result = (logN0 + g * tp) - d * (t - tp);
    }
    
    return result;
    
  }
  
  

  /**
   * Compute mean successful virions
   *
   * @param k hit probability
   * @param V0 inoculation dose
   * @return log mean number of successful virions (based on zero-truncated Poisson distribution with mean k*V0)
   */
  real log_mean_Vs(real k, real V0) {

    real result = log(k*V0) - log_diff_exp(0, -k*V0);
    
    return result;
    
  }

  /**
   * Compute mean successful virions difference
   *
   * @param k hit probability
   * @param V1 inoculation dose one
   * @param V2 inoculation dose two
   * @return difference in log mean number of successful virions (V1-V2)
   */
  real diff_Vs(real k, real V1, real V2) {

    real result = log_mean_Vs(k, V1) - log_mean_Vs(k, V2);

    return result;
    
  }
  
  /**
   * Compute log peak time
   *
   * @param k hit probability
   * @param V0 inoculation dose
   * @param alphat slope parameter for log peak time
   * @param logtp5 intercept parameter for log peak time (log peak time for 10^5 dose)
   * @return log time of peak viral load
   */
  real log_peak_time(real k, real V0, real alphat, real logtp5) {
    
    real result = alphat*diff_Vs(k, V0, 1e+5) + logtp5;

    return result;
    
  }
  
  /**
   * Compute log initial viral load
   *
   * @param k hit probability
   * @param V0 inoculation dose
   * @param alphaN slope parameter for log viral load
   * @param logN5 intercept parameter for log viral load (log initial viral load for 10^5 dose)
   * @return initial log viral load
   */
  real log_init_RNA(real k, real V0, real alphaN, real logN5) {
    
    real result = alphaN*diff_Vs(k, V0, 1e+5) + logN5;
    
    return result;

  }

  /**
   * Compute likelihood, given infected
   *
   * @param log_value data log viral load value
   * @param t time of datapoint
   * @param g exponential growth rate
   * @param d exponential decay rate
   * @param tp time of peak viral load
   * @param logN0 initial log viral load
   * @param sigma viral load measurement standard deviation
   * @param logLOD log lower limit of detection
   * @return likelihood contribution from datapoint with log_value at time t given infected status
   */
  real logL_giv_inf(real log_value, real t,
                    real g, real d,
                    real tp, real logN0,
                    real sigma,
                    real logLOD) {

    // compute model expected log viral load value
    real mean_value = log_viral_load(t, g, d, tp, logN0);
                                       
    
    real result = 0;

    // "-" negative
    if (log_value == -100) {
      // from -infinity to logLOD
      result = normal_lcdf(logLOD | mean_value, sigma);
    }
    
    // "+" positive
    else if (log_value == 100) {
      // from logLOD to +infinity
      result = normal_lccdf(logLOD | mean_value, sigma);
    }
    
    // numeric value
    else {
      result = normal_lpdf(log_value | mean_value, sigma);
    }
    
    return result;

  }
  
  
  /**
   * Compute likelihood, given non-infected
   *
   * @param log_value data log viral load value
   * @param p_false false positive probability
   * @return likelihood contribution from datapoint with log_value given non-infected status
   */
  real logL_giv_uninf(real log_value, real p_false) {
    
    // "-" negative
    if (log_value == -100) {
      return log(1 - p_false);
    }
    
    // "+" positive
    else {
      return log(p_false);
    }
    
  }

}

data {

  int<lower = 0> n_monkeys;
  int<lower = 0> n_datapoints;
  int<lower = 0> n_studies;
  array[n_datapoints] int<lower = 1> individual_ids;
  array[n_datapoints] int<lower = 1> study_ids;
  array[n_monkeys] real<lower = 0> TCID_dose;
  array[n_datapoints] real<lower = 0> day_post_infection;
  array[n_datapoints] real log_viral_load;
  array[n_datapoints] real logLOD;

}

parameters {

  // growth rate: doubling time
  // decay rate: halving time
  real z_mean_log_t_double;
  real z_mean_log_t_half;
  real<lower = 0> z_sigma_log_t_double;
  real<lower = 0> z_sigma_log_t_half;
  // study-specific (hierarchical)
  vector[n_studies] z_t_double;
  vector[n_studies] z_t_half;

  // peak time and log initial RNA
  // for log 10^5 inoculation dose
  real z_mean_log_tp5;
  real z_mean_log_N5;
  real<lower = 0> z_sigma_log_tp5;
  real<lower = 0> z_sigma_log_N5;
  // study-specific (hierarchical)
  vector[n_studies] z_log_tp5;
  vector[n_studies] z_log_N5;
  
  // offsets from log 10^7 inoculation dose
  // used to specify the alpha parameters
  real diff_logtp7_logtp5;
  real z_diff_logN7_logN5;

  // study-specific variation 
  vector[n_studies] z_log_sigma_s;

  // hit probability: ID50
  // false positive probability
  real<lower = (log(log(2)) / (2*log(10))) - (log1m_exp(-10*log(10)) / (2*log(10))) - 2> z_log_ID50;
  real<upper = 4*log(10) / log(1.1)> z_log_p_false;

}

transformed parameters {

  // growth rate: doubling time
  // decay rate: halving time
  real mean_log_t_double = log(0.5) + 0.5*z_mean_log_t_double;
  real mean_log_t_half = log(2) + 0.75*z_mean_log_t_half;
  real<lower = 0> sigma_log_t_double = 0.35*z_sigma_log_t_double;
  real<lower = 0> sigma_log_t_half = 0.35*z_sigma_log_t_half;
  // study-specific (hierarchical)
  vector[n_studies] log_t_double = mean_log_t_double + sigma_log_t_double * z_t_double;
  vector[n_studies] log_t_half = mean_log_t_half + sigma_log_t_half * z_t_half;
  vector[n_studies] g = log(2) ./ exp(log_t_double);
  vector[n_studies] d = log(2) ./ exp(log_t_half);

  // peak time and log initial RNA
  // for log 10^5 inoculation dose
  real mean_log_tp5 = log(2) + 0.5*z_mean_log_tp5;
  real mean_log_N5 = log(1e+3) + log(10)*z_mean_log_N5;
  real<lower = 0> sigma_log_tp5 = 0.35*z_sigma_log_tp5;
  real<lower = 0> sigma_log_N5 = 0.35*z_sigma_log_N5;
  // study-specific (hierarchical)
  vector[n_studies] log_tp5 = mean_log_tp5 + sigma_log_tp5 * z_log_tp5;
  vector[n_studies] log_N5 = mean_log_N5 + sigma_log_N5 * z_log_N5;

  // hit probability: ID50
  // false positive probability
  real<lower = log(log(2)) - log1m_exp(-10*log(10))> log_ID50 = 4*log(10) + log(100)*z_log_ID50;
  real<lower = 1e-10, upper = 1> hit_prob = log(2) / exp(log_ID50) + 1e-10;
  real<upper = 0> log_p_false = -4*log(10) + log(1.1)*z_log_p_false;
  real<lower = 0, upper = 1> p_false = exp(log_p_false);

  // dose-dependence
  // for peak time and log initial RNA
  real diff_logN7_logN5 = log(100) + 1.5*z_diff_logN7_logN5;
  real alpha_t = diff_logtp7_logtp5 / diff_Vs(hit_prob, 1e+7, 1e+5);
  real alpha_N = diff_logN7_logN5 / diff_Vs(hit_prob, 1e+7, 1e+5);

  // study-specific variation
  vector[n_studies] log_sigma_s = log(1) + log(1.5)*z_log_sigma_s;
  vector[n_studies] sigma_s = exp(log_sigma_s);

}

model {

  /////////////////
  // specify priors
  /////////////////

  // growth rate: doubling time
  // decay rate: halving time
  z_mean_log_t_double ~ std_normal();
  z_mean_log_t_half ~ std_normal();
  z_sigma_log_t_double ~ std_normal();
  z_sigma_log_t_half ~ std_normal();
  // study-specific (hierarchical)
  z_t_double ~ std_normal();
  z_t_half ~ std_normal();

  // peak time and log initial RNA
  // for log 10^5 inoculation dose
  z_mean_log_tp5 ~ std_normal();
  z_mean_log_N5 ~ std_normal();
  z_sigma_log_tp5 ~ std_normal();
  z_sigma_log_N5 ~ std_normal();
  // study-specific (hierarchical)
  z_log_tp5 ~ std_normal();
  z_log_N5 ~ std_normal();

  // hit probability: ID50
  // false positive probability
  z_log_ID50 ~ std_normal();
  z_log_p_false ~ std_normal();

  // offsets from log 10^7 inoculation dose
  // used to specify to specify the alpha parameters
  diff_logtp7_logtp5 ~ normal(0, 3);
  z_diff_logN7_logN5 ~ std_normal();

  // study-specific variation
  z_log_sigma_s ~ std_normal();

}
