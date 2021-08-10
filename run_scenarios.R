source("diffused_pharma/experiment_utils.R")

#Some script params
test = TRUE
local = TRUE
n_days = 10
#Some base arguments

#Models to Fit
model_H0 = build_model(dConc ~ (- CL  *  Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(0)},
                       list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_tau"= list("init"= 0.0000001), 
                                    "V"= list("init"= 10, "lower" = 0, "upper" = 100)) 
)

model_H1 = build_model(dConc ~ (- CL  * Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(params$sigma_tau)},
                       list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("V"= list("init"= 10, "lower" = 0, "upper" = 100)) )

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
  Km = 4 # mg per liter
  Vmax = 8
  CL = Vmax / Km # mg per day
  Conc0 =0
  return(list("sigma_eps"=0.001, "sigma_tau"=0, "CL"=CL, "Conc0"=Conc0, "Km"= Km, "V"=V))
}

scenario_base = list(name="base",
                     model_H0=model_H0,
                     model_H1=model_H1,
                     design = design,
                     T_statistic=T_statistic,
                     H0_drift=H0_drift,
                     H1_drift=H1_drift)

if(test)
{
  n_simulations= 10
  n_samples = 50
} else
{
  n_simulations=500
  n_samples=500
  
}
if(local)
{
  path = "C:/Users/roden/Dropbox/Masterarbeit"
} else
{
  path = " "
}
scenarios = list(scenario_base)
for(scenario in scenarios)
{
  run_complete_scenario(scenario, path, n_simulations, n_samples, alpha = 0.05, h=0.02)
}

