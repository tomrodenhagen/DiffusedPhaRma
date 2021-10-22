source("diffused_pharma/experiment_utils.R")
source("diffused_pharma/run_utils.R")
source("scripts/fixtures/models.R")
source("scripts/fixtures/parameters.R")
source("scripts/config.R")
sink(stdout(), type = "message")
options(error = traceback)
options <- commandArgs(trailingOnly = TRUE)
config = load_config(options)
# Some script params
test = config$test
n_days = config$n_days
n_retries = config$n_retries
# Some base arguments
low_noise = 0.01
medium_noise = 0.05
high_noise = 0.1
# Fixtures
T_statistic = function(estimate) {
    return(estimate$sigma_tau)[0]
}
H0_drift = function(t, state, u, params) {
    return(-params$CL * state + u/params$V)
}
H1_drift = function(t, state, u, params) {
    return(-params$CL * params$Km/(params$Km + state) * state + u/params$V)
}
# Models to Fit
model_H0_fixed_low_noise = get_H0_model(low_noise, n_retrys = n_retries)
model_H1_const_diffusion_fixed_low_noise = get_H1_model(low_noise, diffusion_term = "CONSTANT",
    n_retrys = n_retries)
model_H0_fixed_medium_noise = get_H0_model(medium_noise, n_retrys = n_retries)
model_H1_const_diffusion_fixed_medium_noise = get_H1_model(medium_noise, diffusion_term = "CONSTANT",
    n_retrys = n_retries)
model_H0_fixed_high_noise = get_H0_model(high_noise, n_retrys = n_retries)
model_H1_const_diffusion_fixed_high_noise = get_H1_model(high_noise, diffusion_term = "CONSTANT",
    n_retrys = n_retries)
# Dosing and design
get_n_dosis = function(n, d) {
    dosis = list()
    for (i in 0:n) {
        dosis[[paste0("element", i)]] = bolus_dosis(i, d, eps = 0.01)
    }
    return(dosis)
}
get_samples = function(n, kth_dosis = 0) {
    samples = c(0.1, 0.2, 0.3, 0.4, 0.5) + kth_dosis
    for (i in 1:n) {
        samples = c(samples, c(i + 0.1, i - 0.1))
    }
    return(sort(samples))
}
dosis = get_n_dosis(n_days, 250)
design_measure_at_first_dose = list(t_start = 0, t_end = n_days, n_samples = get_samples(n_days),
    dosis = dosis)
design_measure_at_10th_dose = list(t_start = 0, t_end = n_days, n_samples = get_samples(n_days,
    10), dosis = dosis)
models = list(list(H0 = model_H0_fixed_low_noise, H1 = model_H1_const_diffusion_fixed_low_noise,
    desc = "Constant diffusion term with fixed low sigma", shortcut = "const_diff_fixed_low_noise"),
    list(H0 = model_H0_fixed_high_noise, H1 = model_H1_const_diffusion_fixed_high_noise,
        desc = "Constant diffusion term with fixed high sigma", shortcut = "const_diff_fixed_high_noise"),
    list(H0 = model_H0_fixed_medium_noise, H1 = model_H1_const_diffusion_fixed_medium_noise,
        desc = "Constant diffusion term with fixed medium sigma", shortcut = "const_diff_fixed_medium_noise"))
parameter_samplings = list(list(sampling = get_parameter_sampling(low_noise), desc = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is low variance",
    shortcut = "low_noise"), list(sampling = get_parameter_sampling(medium_noise),
    desc = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is medium variance",
    shortcut = "medium_noise"), list(sampling = get_parameter_sampling(high_noise),
    desc = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is high variance",
    shortcut = "high_noise"))
designs = list(list(design = design_measure_at_first_dose, desc = "Dosing at every day and measure shortly after and shortly before.
               Additionaly measure some more points after the first dose",
    shortcut = "measure_first_cycle"))
config$name = "fixed_noise"
run_scenarios(models, designs, parameter_samplings, config)
