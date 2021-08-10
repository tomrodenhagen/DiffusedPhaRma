source("diffused_pharma/fit.R")
source("diffused_pharma/sim.R")
library(parallel)

run_test = function(data, model_H0, model_H1, T_statistic, n_simulations, design, alpha=0.95,h=0.1)
{ 
  estimated_params_H0 = model_H0$estimate(data)

  estimated_params_H1 = model_H1$estimate(data)
  T_samples = c()
  if(.Platform$OS.type == "unix")
  {
    numCores = detectCores()
  }
  else
  {
    numCores <- 1
  }
  generate_T_sample = function(dummy)
  {   
      data_sim = model_H0$simulate(estimated_params_H0, design$t_start, design$t_end,design$n_samples, h, design$dosis)
      plot(data_sim[["t"]], data_sim[["ConcObserved"]])
      estimate_sim = model_H1$estimate(data_sim)
      
      return(T_statistic(estimate_sim))		
  }
  T_samples <-mclapply(1:n_simulations, generate_T_sample, mc.cores = numCores)
  T_samples = unlist(T_samples, use.names=FALSE)
  emp_quantile = quantile(T_samples, 1 - alpha)
  print(emp_quantile)
  print(T_statistic(estimated_params_H1))
  return(list(rejected = emp_quantile < T_statistic(estimated_params_H1)))
}
visualize_setting= function(drift, diffusion, model_H0, model_H1, T_statistic, sample_params, design,h)
{
  sampled_params = sample_params()
  data_obs = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end, design$n_samples, h=h, design$dosis)
  data_unobs = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end,
                              1000, h=h, design$dosis)
  
  plot(data_obs[["t"]], data_obs[["ConcObserved"]])
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green")
  estimated_params_H0 = model_H0$estimate(data_obs)
  estimated_params_H1 = model_H1$estimate(data_obs)
  print(estimated_params_H0)
  print(estimated_params_H1)
  data_H0 = model_H0$simulate(estimated_params_H0, design$t_start, design$t_end,design$n_samples,  h=h, design$dosis)
  data_H1 = model_H1$simulate(estimated_params_H1, design$t_start, design$t_end,design$n_samples,  h=h, design$dosis)
  lines(data_H0[["t"]], data_H0[["ConcObserved"]], pch=2, col="red")
  lines(data_H1[["t"]], data_H1[["ConcObserved"]], pch=3, col="blue")
  
}
library(progress)
run_simulation_study = function(drift, diffusion, model_H0, model_H1, T_statistic, sample_params, design, n_param_samples, n_simulations, h=0.1)
{ 
  visualize_setting(drift, diffusion, model_H0, model_H1, T_statistic, sample_params, design,h)
  rec <- as.data.frame(matrix(0, ncol = length(sample_params()), nrow = n_param_samples))
  res_vec = c()
  names(rec) = names(sample_params())
  pb = progress_bar$new(total = n_param_samples)
  for(i in 1:n_param_samples)
  { pb$tick()
    sampled_params = sample_params()
   
    rec[i,] = unlist(sampled_params, use.names = FALSE) 
    data = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end, design$n_samples, h=h, design$dosis)
    
    #Data should look like real world data
    res = run_test(data, model_H0, model_H1, T_statistic, n_simulations, design,h)
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


