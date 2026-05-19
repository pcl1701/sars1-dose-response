# functions for SARS-CoV-1 dose-response analysis

# unit conversions -------------------------------------------------------------

# log10 to log
log10_to_log <- function (log10_value) {
  return (log10_value / log10(exp(1)))
}

# log to log10
log_to_log10 <- function (log_value) {
  return (log_value / log(10))
}

# dose-response ----------------------------------------------------------------

# dose-response probability of infection
# INPUTS:
# (hit_prob): hit probability
# (TCID_dose): TCID50 inoculation dose
# OUTPUTS:
# returns probability of infection given (TCID_dose) dose
DR_prob <- function(hit_prob, TCID_dose) {
  return (1 - exp(-hit_prob * TCID_dose))
}

# viral kinetics ---------------------------------------------------------------

# exponential growth and decay
# INPUTS:
# (t): time (days)
# (g_i): growth rate
# (d_i): decay rate
# (t_p): peak time (days)
# (logN_0): initial log viral load
# OUTPUTS:
# returns log viral load value at time (t)
exp_growth_decay <- function (t,
                              g_i, d_i,
                              t_p, logN_0) {
  # growth phase
  if (t < t_p) {
    return (logN_0 + g_i * t)
  }
  # decay phase
  else {
    return ((logN_0 + g_i * t_p) - d_i * (t - t_p))
  }
}

# successful virions
# INPUTS:
# (k): hit probability
# (V_0): initial inoculation dose, units of TCID50
# OUTPUTS:
# returns log mean number of successful virions, units of log TCID50
# given there was at least one successful virion
# (log mean of zero-truncated Poisson distribution with mean k*V_0)
log_mean_V_s <- function(k, V_0) {
  return (log(k*V_0) - log(exp(0) - exp(-k*V_0)))
}

# dose-dependence, log peak time
# INPUTS:
# (k): hit probability
# (V_0): TCID50 inoculation dose
# (alpha_t): slope parameter for log peak time
# (beta_t): intercept parameter for log peak time
# OUTPUTS:
# returns log peak time, units of log days
peak_time <- function(k, V_0, alpha_t, beta_t) {
  return (exp(alpha_t*(log_mean_V_s(k, V_0) - log_mean_V_s(k, 1e+5)) + beta_t))
}

# dose-dependence, URT log initial viral load
# INPUTS:
# (k): hit probability
# (V_0): TCID50 inoculation dose
# (alpha_N): slope parameter for log viral load
# (beta_N): intercept parameter for log viral load
# OUTPUTS:
# returns URT log initial viral load
log_init_RNA <- function(k, V_0, alpha_N, beta_N) {
  return (alpha_N*(log_mean_V_s(k, V_0) - log_mean_V_s(k, 1e+5)) + beta_N);
}
