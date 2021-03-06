source("diffused_pharma/experiment_utils.R")
source("diffused_pharma/run_utils.R")
source("scripts/fixtures/models.R")
source("scripts/fixtures/parameters.R")
source("scripts/fixtures/design.R")
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
model_H0 = get_H0_model(n_retrys = n_retries)
model_H1_const_diffusion = get_H1_model(diffusion_term = "CONSTANT", n_retrys = n_retries)
model_H1_linear_diffusion = get_H1_model(diffusion_term = "LINEAR", n_retrys = n_retries)
dosis = get_n_dosis(n_days, 250)
design_measure_at_first_dose_many_samples = list(t_start = 0, t_end = n_days, n_samples = get_samples(n_days,
    sampling = "MANY"), dosis = dosis)
design_measure_at_first_dose_few_samples = list(t_start = 0, t_end = n_days, n_samples = get_samples(n_days,
    sampling = "FEW"), dosis = dosis)
models = list(list(H0 = model_H0, H1 = model_H1_const_diffusion, desc = "Constant diffusion term with sigma epsilon to fit",
    shortcut = "const_diff"), list(H0 = model_H0, H1 = model_H1_linear_diffusion,
    desc = "Linear diffusion term with sigma epsilon to fit", shortcut = "linear_diff"))
parameter_samplings = list(list(sampling = get_parameter_sampling(low_noise), desc = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is low variance",
    shortcut = "low_noise"), list(sampling = get_parameter_sampling(medium_noise),
    desc = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is medium variance",
    shortcut = "medium_noise"), list(sampling = get_parameter_sampling(high_noise),
    desc = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is high variance",
    shortcut = "high_noise"))
designs = list(list(design = design_measure_at_first_dose_many_samples, desc = "Dosing at every day and measure shortly after and shortly before.
               Additionaly measure some many points after the first dose",
    shortcut = "measure_first_cycle_many_samples"), list(design = design_measure_at_first_dose_few_samples,
    desc = "Dosing at every day and measure shortly after and shortly before.
\t\t            Additionaly measure some few points after the first dose",
    shortcut = "measure_first_cycle_few_samples"))
config$name = "diff_measurements"
run_scenarios(models, designs, parameter_samplings, config)
