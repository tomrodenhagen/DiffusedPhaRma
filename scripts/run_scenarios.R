source("diffused_pharma/experiment_utils.R")
sink(stdout(), type="message")
options(error=traceback)
#Some script params
test = TRUE
local = FALSE
n_days = 15
#Some base arguments

#Models to Fit
sigma_eps_fixed = 0.05
model_H0 = build_model(dConc ~ (- CL  *  Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(0)},
                       params=list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_eps"= list("init"= sigma_eps_fixed),
				    "sigma_tau"= list("init"= 0.0000001),
                                   "V"= list("init"= 10, "lower" = 0, "upper" = 100)) 
)
model_H1_const_diff = build_model(dConc ~ (- CL  *  Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(0)},
                       params=list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_eps"= list("init"= sigma_eps_fixed),
                                   "V"= list("init"= 10, "lower" = 0, "upper" = 100)) 
)

model_H1_linear_diff = build_model(dConc ~ (- CL - sigma_tau**2 * 0.5  +  D / V / (exp(Conc) )  ) * dt + sigma_tau * dw1, 
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V  )},
                       diffusion = function(t, state, u , params){return(params$sigma_tau * state  )},
                       params=list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_eps"= list("init"= sigma_eps_fixed),
				    "V"= list("init"= 10, "lower" = 0, "upper" = 100),
				    "CL"= list("init"= 1, "lower" = 0.1, "upper" = 10),
				    "Conc0"= list("init"= -5),
				    "sigma_tau"= list("init"= 0.2, "lower" = 0.01, "upper" = 1)
				   ),
		       transf = function(x){return(exp(x))},
                       observation_equation = ConcObserved ~ exp(Conc)
                       
)
 






T_statistic = function(estimate){return(estimate$sigma_tau)[0]}
#Models to simulate
H0_drift =function(t, state, u, params){return(-params$CL * state + u / params$V) }
H1_drift =function(t, state, u, params){return(-params$CL * params$Km / (params$Km + state) * state + u / params$V) }


#Dosing and design

get_n_dosis = function(n, d)
{ dosis = list()
for(i in 0:n)
{
  dosis[[paste0("element", i)]] = bolus_dosis(i , d, eps=0.01)
}
return(dosis)
}

get_samples = function(n)
{ samples = c(0.1, 0.2, 0.3, 0.4, 0.5)
for(i in 1:n)
{
  samples = c(samples, c(i + 0.1, i - 0.1 ))
}
return(samples)
}
dosis = get_n_dosis(n_days, 250)
design = list(t_start=0, t_end = n_days, n_samples=get_samples(n_days), dosis = dosis )

#Choice of parameters
sample_params = function()
{ V = 50
  Km = sample(c(2,3,4,5,6,7),1)  # mg per liter
  Vmax = 7
  CL = Vmax / 4 # mg per day
  Conc0 = 0.03
  return(list("sigma_eps"=sigma_eps_fixed, "sigma_tau"=0, "CL"=CL, "Conc0"=Conc0, "Km"= Km, "V"=V))
}

scenario_1 = list(name="const_diff_test",
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




if(test)
{
  n_simulations_typ1=200
  n_simulations_typ2=200
  n_samples=100

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

