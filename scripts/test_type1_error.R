source("diffused_pharma/experiment_utils.R")

model_H0 = build_model(dConc ~ (- CL  * Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(0)},
                       list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_tau"= list("init"= 0.0000001))
)

model_H1 = build_model(dConc ~ (- CL  * Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(params$sigma_tau)},
                       list("Conc0","CL", "sigma_eps","sigma_tau", "V"))

T_statistic = function(estimate){return(estimate$sigma_tau)[0]}

design = list(t_start=0, t_end = 10, n_samples=seq(0.1,10,length.out=20), dosis = list(bolus_dosis(0.05, 2), 
                                                                                       bolus_dosis(2.05, 2), 
                                                                                       bolus_dosis(4.05, 2), bolus_dosis(6.05, 2)) )

prior_H0 = function()
{ 
  CL = 1
  Conc0 =0
  return(list("sigma_eps"=0.001, "sigma_tau"=0, "CL"=CL, "Conc0"=Conc0, "Km"= 2, "V"=2))
  }

rec = run_simulation_study(function(t, state, u, params){return(-params$CL * state + u / params$V)},
                     diffusion = function(t, state, u , params){return(0)},
                     model_H0, model_H1, T_statistic, prior_H0, design, 50,100 )

eval_simulation(rec, "test")


