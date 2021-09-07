source("diffused_pharma/experiment_utils.R")
source("scripts/fixtures/models.R")
sink(stdout(), type="message")
options(error=traceback)
#Some script params
test = TRUE
local = .Platform$OS.type != "unix"
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
return(sample_params)
}


models = list( list("H0" = model_H0, "H1"= model_H1_const_diffusion, 
	      "desc" = "Constant diffusion term with sigma epsilon to fit", "shortcut" = "const_diff") ,
	      list("H0" = model_H0, "H1"= model_H1_linear_diffusion,
	      "desc" = "Linear diffusion term with sigma epsilon to fit", "shortcut" = "linear_diff"))

parameter_samplings = list(list("sampling" = get_parameter_sampling(low_noise),
                          "desc" = "Uniform sampling of Km between 3 and 9,  CL=1.75, V=50 and sigma_eps is low_var",
			                    "shortcut" = "low_noise"))
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
scenarios = list()
for( m in models)
{
  for(s in parameter_samplings)
  {
    for (d in designs)
    { 
      name = paste(m$shortcut, s$shortcut, d$shortcut, suffix, sep="+")
      desc = paste(m$desc, s$desc, d$desc,  sep="::")
      scenario = list(name=name,
                      description=desc,
                      model_H0=m$H0,
                      model_H1=m$H1,
                      design = d$design,
                      sample_params=s$sampling,
                      T_statistic=T_statistic,
                      H0_drift=H0_drift,
                      H1_drift=H1_drift)
      scenarios[[name]] = scenario
    }
  }
  
}

if(test)
{
  n_simulations_typ1=50
  n_simulations_typ2=50
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

for(scenario in scenarios)
{ print(scenario$name)
  res = run_complete_scenario(scenario, path, n_simulations_typ1, n_simulations_typ2, n_samples, alpha = 0.05, h=0.02, TRUE)
 
  print(eval_simulation(res[["type1"]]))
  print(eval_simulation(res[["type2"]]))
}

