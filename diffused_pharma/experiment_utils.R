source("diffused_pharma/fit.R")
source("diffused_pharma/sim.R")
library(parallel)
library(jsonlite)
library(unixtools)
run_test = function(data, model_H0, model_H1, T_statistic, n_simulations, design, alpha=0.05,h=0.1)
{ 
  estimated_params_H0 = model_H0$estimate(data)

  estimated_params_H1 = model_H1$estimate(data)
  T_samples = c()
  generate_T_sample = function(dummy)
  {   
      data_sim = model_H0$simulate(estimated_params_H0, design$t_start, design$t_end,design$n_samples, h, design$dosis)
      estimate_sim = model_H1$estimate(data_sim)     
      return(T_statistic(estimate_sim))		
  }
  T_samples <-mclapply(1:n_simulations, generate_T_sample, mc.cores = 1)
  T_samples = unlist(T_samples, use.names=FALSE)
  
  emp_quantile = quantile(T_samples, 1 - alpha)
  return(list(rejected = emp_quantile < T_statistic(estimated_params_H1)))
}
visualize_setting= function(drift, diffusion, model_H0, model_H1, T_statistic, sample_params, design, h, path, name, draw_vmax=FALSE)
{
  sampled_params = sample_params()
  data_obs = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end, design$n_samples, h=h, design$dosis)
  data_unobs = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end,
                              1000, h=h, design$dosis)
  estimated_params_H0 = model_H0$estimate(data_obs)
  estimated_params_H1 = model_H1$estimate(data_obs)
  pathH0 = file.path(path, paste(name,"visH0.png",sep="_"))
  jpeg(file=pathH0)
  plot(data_obs[["t"]], data_obs[["ConcObserved"]],
       main="Sample experiment",
       ylab="Drug concentration in mg /l",
       xlab="time in days", pch=19,col="blue",cex=1.5) 
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green")
  
  
  data_H0 = model_H0$simulate(estimated_params_H0, design$t_start, design$t_end,300,  h=h, design$dosis)
  lines(data_H0[["t"]], data_H0[["ConcObserved"]], pch=2, col="red")
  legend( x= "bottomright", legend=c("Observed points", "True", "Predicted from Model H0"),
  col=c("blue","green", "red"), lty=1:2, cex=1,pch=19)
  if(draw_vmax)
    {
    abline(h=sampled_params$Vmax, col="brown")
    }
  dev.off()
  
  pathH1 = file.path(path, paste(name,"visH1.png",sep="_"))
  jpeg(file=pathH1)
  plot(data_obs[["t"]], data_obs[["ConcObserved"]],
       main="Sample experiment",
       ylab="Drug concentration in mg /l",
       xlab="time in days", pch=19,col="blue",cex=1.5) 
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green")
 
 
  data_H1 = model_H1$simulate(estimated_params_H1, design$t_start, design$t_end,300,  h=h, design$dosis)
  lines(data_H1[["t"]], data_H1[["ConcObserved"]], pch=3, col="blue")
  if(draw_vmax)
  {
    abline(h=sampled_params$Vmax, col="brown")
  }
  legend( x= "bottomright", legend=c("Observed points", "True", "Predicted from Model H1"),
          col=c("blue","green", "blue"), lty=1:2, cex=1)
  
  dev.off()
  
}
library(progress)
run_simulation_study = function(drift, diffusion, model_H0, model_H1, T_statistic, sample_params, design, n_param_samples, n_simulations, h=0.1)
{ 
  
  rec <- as.data.frame(matrix(0, ncol = length(sample_params()), nrow = n_param_samples))
  
  names(rec) = names(sample_params())
  res_vec = c()  
  pb = progress_bar$new(total = n_param_samples)
  if(.Platform$OS.type == "unix")
  {
    numCores = detectCores()
  }
  else
  {
    numCores <- 1
  } 
  run_test_wrapper = function(dummy)
  { tmp_path = file.path("/tmp", paste("tmp", dummy, sep="_") )
    
    dir.create(tmp_path)
    set.tempdir(tmp_path)
    sampled_params = sample_params() 
    data = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end, design$n_samples, h=h, design$dosis)
    #Data should look like real world data
    model_H0_copy=copy_model(model_H0)
    model_H1_copy=copy_model(model_H1)
    res = run_test(data, model_H0_copy, model_H1_copy, T_statistic, n_simulations, design,h)
    
    return(res)
  }
  
  res_vec = mclapply(1:n_param_samples, run_test_wrapper, mc.cores = numCores)
  res = unlist(res_vec, use.names=FALSE)
  print(res)
  rec[["test_rejected"]] = res
  
  return(rec)
 
}
eval_simulation = function(rec, path=NULL, name=NULL )
{ 
  res_evaluated = ("Percentages of rejected tests"=mean(rec$test_rejected) ) 
  if(!is.null(path))
  {save(rec, file= file.path(path, paste(name , "complete_res.csv", sep = "_") ) )
   res_evaluated <- toJSON(res_evaluated)
   write(res_evaluated, file.path(path,paste(name , "eval.csv", sep = "_")))
  }
  return(res_evaluated)
 }
run_complete_scenario = function(scenario, path = "C:/Users/roden/Dropbox/Masterarbeit", n_simulations = 500, n_samples=500, alpha = 0.05, h=0.02, fresh=TRUE)
{ scenario_folder = file.path(path, scenario$name)
  if(dir.exists(scenario_folder))
  { if(fresh)
    {
      unlink(scenario_folder,rec=TRUE)
      dir.create(scenario_folder)
    }
  } else
  {
    dir.create(scenario_folder)
  }

  
  #some standart parameters, we might change that at a later point
  diffusion = function(t, state, u , params){return(0)}
  #Visualize 
  visualize_setting(scenario$H0_drift,
                    diffusion, 
                    scenario$model_H0, 
                    scenario$model_H1,
                    scenario$T_statistic,
                    scenario$sample_params,
                    scenario$design,
                    h, 
                    scenario_folder, 
                    "vis_type1") 
  visualize_setting(scenario$H1_drift,
                    diffusion, 
                    scenario$model_H0, 
                    scenario$model_H1,
                    scenario$T_statistic,
                    scenario$sample_params,
                    scenario$design,
                    h, 
                    scenario_folder,
                    "vis_type2",
                    draw_vmax = TRUE) 
  #Type 1 Error
  res_type_1 = run_simulation_study(scenario$H0_drift,
                             diffusion,
                             scenario$model_H0, 
                             scenario$model_H1, 
                             scenario$T_statistic, 
                             scenario$sample_params, 
                             scenario$design, n_simulations, n_samples, alpha)
  eval_simulation(res_type_1, scenario_folder, "typ1_res")
  #Type 2 error
  res_type_2 = run_simulation_study(scenario$H1_drift,
                                    diffusion,
                                    scenario$model_H0, 
                                    scenario$model_H1, 
                                    scenario$T_statistic, 
                                    scenario$sample_params, 
                                    scenario$design, n_simulations, n_samples, alpha)
  eval_simulation(res_type_2, scenario_folder, "typ2_res")
  return(list(type1=res_type_1, type2=res_type_2))
}
  


