source("diffused_pharma/experiment_utils.R")
source("diffused_pharma/run_utils.R")
source("scripts/fixtures/models.R")
sink(stdout(), type="message")
options(error=traceback)

#Some script params
test = TRUE
local = .Platform$OS.type != "unix"
n_days = 15
if(test)
{
 n_retries=4
} else
{
 n_retries=10
}
#Some base arguments

low_noise = 0.01
high_noise = 0.05

#Fixtures
T_statistic = function(estimate){return(estimate$sigma_tau)[0]}

H0_drift =function(t, state, u, params){return(-params$CL * state + u / params$V) }
H1_drift =function(t, state, u, params){return(-params$CL * params$Km / (params$Km + state) * state + u / params$V) }



#Models to Fit
model_H0 =  get_H0_model(n_retrys=n_retries)
model_H1_const_diffusion = get_H1_model(diffusion_term="CONSTANT", n_retrys=n_retries)
model_H1_linear_diffusion = get_H1_model(diffusion_term="LINEAR", n_retrys=n_retries)
 
model_H0_fixed_low_noise =  get_H0_model(low_noise, n_retrys=n_retries)
model_H1_const_diffusion_fixed_low_noise = get_H1_model(low_noise, diffusion_term="CONSTANT", n_retrys=n_retries)
model_H1_linear_diffusion_fixed_low_noise = get_H1_model(low_noise, diffusion_term="LINEAR", n_retrys=n_retries)
 
model_H0_fixed_high_noise =  get_H0_model(high_noise, n_retrys=n_retries)
model_H1_const_diffusion_fixed_high_noise = get_H1_model(high_noise, diffusion_term="CONSTANT", n_retrys=n_retries)
model_H1_linear_diffusion_fixed_high_noise = get_H1_model(high_noise, diffusion_term="LINEAR", n_retrys=n_retries)
 



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
return(sample_params)
}


models = list( list("H0" = model_H0, "H1"= model_H1_const_diffusion, 
	       "desc" = "Constant diffusion term with sigma epsilon to fit", "shortcut" = "const_diff") ,
	       list("H0" = model_H0, "H1"= model_H1_linear_diffusion,
	       "desc" = "Linear diffusion term with sigma epsilon to fit", "shortcut" = "linear_diff"),
	       list("H0" = model_H0_fixed_low_noise, "H1"= model_H1_const_diffusion_fixed_low_noise,
		   "desc" = "Constant diffusion term with fixed low sigma", "shortcut" = "const_diff_fixed_low_noise"),
	       list("H0" = model_H0_fixed_low_noise, "H1"= model_H1_linear_diffusion_fixed_low_noise,
		   "desc" = "Linear diffusion term with fixed low sigma", "shortcut" = "linear_diff_fixed_low_noise"),
	      list("H0" = model_H0_fixed_high_noise, "H1"= model_H1_const_diffusion_fixed_high_noise,
                   "desc" = "Constant diffusion term with fixed high sigma", "shortcut" = "const_diff_fixed_high_noise"),
	       list("H0" = model_H0_fixed_high_noise, "H1"= model_H1_linear_diffusion_fixed_high_noise,
                   "desc" = "Linear diffusion term with fixed high sigma", "shortcut" = "linear_diff_fixed_high_noise"))

parameter_samplings = list( list("sampling" = get_parameter_sampling(low_noise),
                                 "desc" = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is high variance",
			         "shortcut" = "low_noise"),
			    list("sampling" = get_parameter_sampling(high_noise),
                                 "desc" = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is high variance",
                                 "shortcut" = "high_noise")
			   )
designs = list( list("design"=design_measure_at_first_dose,
               "desc"= "Dosing at every day and measure shortly after and shortly before.
               Additionaly measure some more points after the first dose",
               "shortcut" = "measure_first_cycle"))

if(test)
{
  suffix = "test"
} else
{
  suffix="full"
}

if(test)
{
  n_simulations_typ1= 200
  n_simulations_typ2=200
  n_samples=200

} else
{
  n_simulations_typ1=400
  n_simulations_typ2=5000
  n_samples= 300
  
}
if(local)
{
  base_path = "C:/Users/roden/Documents/data"
} else
{
  base_path = "/localhome/tr/runs"
}
config = list(name = "generic_combination", 
          fresh = TRUE,
	  base_path = base_path,
	  n_simulations_typ1 = n_simulations_typ1,
	  n_simulations_typ2=n_simulations_typ2,
	  n_samples=n_samples,
          alpha = 0.05,
          h = 0.02)  
run_scenarios(models, designs, parameter_samplings, config)
