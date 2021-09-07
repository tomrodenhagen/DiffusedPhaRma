source("diffused_pharma/experiment_utils.R")
source("scripts/fixtures/models.R")
sink(stdout(), type="message")
options(error=traceback)
#Some script params
test = TRUE
local = FALSE
n_days = 15
#Some base arguments

low_noise = 0.01
high_noise = 0.05

#Fixtures
T_statistic = function(estimate){return(estimate$sigma_tau)[0]}

H0_drift =function(t, state, u, params){return(-params$CL * state + u / params$V) }
H1_drift =function(t, state, u, params){return(-params$CL * params$Km / (params$Km + state) * state + u / params$V) }



#Models to Fit
model_H0 =  get_H0_model()
model_H1_const_diffusion = get_H1_model(diffusion_term="CONST")
model_H1_linear_diffusion = get_H1_model(diffusion_term="LINEAR")
 
model_H0_fixed_low_noise =  get_H0_model(low_noise)
model_H1_const_diffusion_fixed_low_noise = get_H1_model(low_noise, diffusion_term="CONST")
model_H1_linear_diffusion_fixed_low_noise = get_H1_model(low_noise, diffusion_term="LINEAR")
 
model_H0_fixed_high_noise =  get_H0_model(high_noise)
model_H1_const_diffusion_high_noise = get_H1_model(high_noise, diffusion_term="CONST")
model_H1_linear_diffusion_high_noise = get_H1_model(high_noise, diffusion_term="LINEAR")
 



#Dosing and design

get_n_dosis = function(n, d)
{ dosis = list()
  for(i in 0:n)
  {
    dosis[[paste0("element", i)]] = bolus_dosis(i , d, eps=0.01)
  }
  return(dosis)
}

get_samples = function(n, kth_dosis=0)
{ samples = c(0.1, 0.2, 0.3, 0.4, 0.5) + kth_dosis
for(i in 1:n)
{
  samples = c(samples, c(i + 0.1, i - 0.1 ))
}
return(sort(samples))
}
dosis = get_n_dosis(n_days, 250)
design_measure_at_first_dose = list(t_start=0, t_end = n_days, n_samples=get_samples(n_days), dosis = dosis )
design_measure_at_10th_dose = list(t_start=0, t_end = n_days, n_samples=get_samples(n_days, 10), dosis = dosis )

#Choice of parameters
get_parameter_sampling = function(sigma_eps)
{
sample_params = function()
{ V = 50
  Km = sample(c(3,4,5,6,7,8,9),1)  # mg per liter
  Vmax = 7
  CL = Vmax / 4 # mg per day
  Conc0 = 0
  return(list("sigma_eps"=sigma_eps, "sigma_tau"=0, "CL"=CL, "Conc0"=Conc0, "Km"= Km, "V"=V))
}
return sample_params
}
scenario_1 = list(name="constant_diffusion",
		     description="Constant diffusionterm. Fixed sigma_eps=0.05 in the model. 
		     		  Choosed Km between 2 and 7",
                     model_H0=model_H0,
                     model_H1=model_H1_const_diff,
                     design = design,
                     sample_params=sample_params,
                     T_statistic=T_statistic,
                     H0_drift=H0_drift,
                     H1_drift=H1_drift)

scenario_2 = list(name="linear_diff_test",
		     description="Linear diffusionterm. Fixed sigma_eps=0.05 in the model. 
		     		  Choosed Km between 2 and 7",
                     model_H0=model_H0,
                     model_H1=model_H1_linear_diff,
                     design = design,
                     sample_params=sample_params,
                     T_statistic=T_statistic,
                     H0_drift=H0_drift,
                     H1_drift=H1_drift)

models = c( c("H0" = model_H0, "H1"= model_H1_constant_diffusion, 
	      "desc" = "Constant diffusion term with sigma epsilon to fit", "shortcut" = "const_diff") ,
	    c("H0" = model_H0, "H1"= model_H1_constant_diffusion,
	      "desc" = "Linear diffusion term with sigma epsilon to fit", "shortcut" = "linear_diff"))

parameter_samplings = c(c("sampling" = get_parameter_sampling(low_var), "desc" = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is low_var",
			  "shortcut" = 


if(test)
{
  n_simulations_typ1=100
  n_simulations_typ2=200
  n_samples=50

} else
{
  n_simulations_typ1=400
  n_simulations_typ2=5000
  n_samples= 300
  
}
if(local)
{
  path = "C:/Users/roden/Documents/data"
} else
{
  path = "/localhome/tr"
}
scenarios = list(scenario_1, scenario_2)
for(scenario in scenarios)
{
  res = run_complete_scenario(scenario, path, n_simulations_typ1, n_simulations_typ2, n_samples, alpha = 0.05, h=0.02, TRUE)
  print(scenario$name)
  print(eval_simulation(res[["type1"]]))
  print(eval_simulation(res[["type2"]]))
}

