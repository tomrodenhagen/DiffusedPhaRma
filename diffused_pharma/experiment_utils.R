source("diffused_pharma/fit.R")
source("diffused_pharma/sim.R")
library(parallel)

run_test = function(data, model_H0, model_H1, T_statistic, n_simulations, alpha=0.95)
{
  estimated_params_H0 = model_H0$estimate(data)
  estimated_params_H1 = model_H1$estimate(data)
  T_samples = c()
  numCores <- detectCores()
  generate_T_sample = function(dummy)
  {   
      data_sim = model_H0$simulate(estimated_params_H0, attr(data,"t_start"), attr(data,"t_end"), attr(data,"n_samples"), h=0.1, attr(data,"dosis"))
      estimate_sim = model_H1$estimate(data_sim)
      return(T_statistic(estimate_sim))		
  }
  T_samples <-mclapply(1:n_simulations, generate_T_sample, mc.cores = numCores)
  T_samples = unlist(T_samples, use.names=FALSE)
  emp_quantile = quantile(T_samples, alpha)
  return(list(rejected = emp_quantile < T_statistic(estimated_params_H1)))
}
library(progress)
run_simulation_study = function(drift, diffusion, model_H0, model_H1, T_statistic, sample_params, design, n_param_samples, n_simulations)
{ 
  
  rec <- as.data.frame(matrix(0, ncol = length(sample_params()), nrow = n_param_samples))
  res_vec = c()
  names(rec) = names(sample_params())
  pb = progress_bar$new(total = n_param_samples)
  for(i in 1:n_param_samples)
  { pb$tick()
    sampled_params = sample_params()
   
    rec[i,] = unlist(sampled_params, use.names = FALSE) 
    data = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end, design$n_samples, h=0.1, design$dosis)
    res = run_test(data, model_H0, model_H1, T_statistic, n_simulations)
    res_vec = c(res_vec, as.integer(res$rejected) )
  }
  rec[["test_rejected"]] = res_vec
  return(rec)
 
}
eval_simulation = function(rec, name )
{ 
  print(sprintf("Percentages of rejected tests: %s", mean(rec$test_rejected) ) )
  for(param in colnames(rec))
  {
    if(param=="sigma_eps" | param=="test_rejected")
    {
      next
    }
    #plot(ksmooth(rec[[param]], rec[["test_rejected"]]),title=param,xlab=param,ylab="test_rej")
  }
}


