source("diffused_pharma/experiment_utils.R")

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


sample_params = function()
{ V = 50
  Km = 1000 # mg per liter
  Vmax = 100
  CL = Vmax / Km # mg per day
  Conc0 =0
  return(list("sigma_eps"=0.001, "sigma_tau"=0, "CL"=CL, "Conc0"=Conc0, "Km"= Km, "V"=V))
}


get_n_dosis = function(n, d)
{ dosis = list()
  for(i in 0:n)
  {
    dosis[[paste0("element", i)]] = bolus_dosis(i , d, eps=0.01)
  }
  return(dosis)
}

get_samples = function(n)
{ samples = c(0.1)
  for(i in 1:n)
  {
    samples = c(samples, c(i + 0.1, i - 0.1 ))
  }
  return(samples)
}

n_days = 15


dosis = get_n_dosis(n_days, 250)

design = list(t_start=0, t_end = n_days, n_samples=get_samples(n_days), dosis = dosis )
#visualize_setting(function(t, state, u, params){return(-params$CL * params$Km / (params$Km + state) * state + u / params$V)},
#                                                          diffusion = function(t, state, u , params){return(0)},
#                                                       model_H0, model_H1, T_statistic, sample_params = sample_params, design,  0.05 )              

H0_drift =function(t, state, u, params){return(-params$CL * state + u / params$V) }
H1_drift =function(t, state, u, params){return(-params$CL * params$Km / (params$Km + state) * state + u / params$V) }

rec = run_simulation_study(H0_drift,
                                         diffusion = function(t, state, u , params){return(0)},
                                         model_H0, model_H1, T_statistic, sample_params, design, 10,100, 0.05 )
              
eval_simulation(rec, "test")