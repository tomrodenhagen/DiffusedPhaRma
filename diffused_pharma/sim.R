

next_step_sde = function(drift, diffusion, t, state, u ,params, h)
{ #Euler schema
  drift_value = drift(t, state, u, params)
  diffusion_value = diffusion(t, state,u, params)
  state + drift_value * h + diffusion_value * h**0.5 * rnorm(1,0,1)
}
library(tidyr)
get_t_grid = function(t_start, t_end, n_samples, h, list_of_dosis)
{ 
  if(length(n_samples) > 1)
  { 
    output_samples = n_samples
    
  }
  else
  { 
    output_samples = seq(t_start, t_end,length.out=n_samples)
  }
  simulation_samples = seq(t_start, t_end, h)
  dosis_samples = c()
  rate = c()
  
  for(dosis in list_of_dosis)
  {
    dosis_samples = c(dosis_samples, dosis$start)
    dosis_samples = c(dosis_samples, dosis$end)
    rate = c(rate, dosis$rate)
    rate = c(rate, 0)
  }
  
  output_samples_df = data.frame("t" = output_samples, 
                                 "output" = rep(TRUE, length(output_samples)),
                                 "D" = rep(NA, length(output_samples)),
                                 "D_sample" = rep(FALSE, length(output_samples))               
                                 )
  simulation_samples_df = data.frame("t" = simulation_samples, 
                                 "output" = rep(FALSE, length(simulation_samples)),
                                 "D" = rep(NA, length(simulation_samples)),
                                 "D_sample" = rep(FALSE, length(simulation_samples))
  )
  
  dosis_samples_df = data.frame("t" = dosis_samples, 
                                     "output" = rep(FALSE, length(dosis_samples)),
                                     "D" = rate,
                                     "D_sample" = rep(TRUE, length(dosis_samples))
                                )
  combined_df = rbind(output_samples_df, simulation_samples_df, dosis_samples_df)
  combined_df = combined_df[order(combined_df$t),]
  combined_df = fill(combined_df, D)
  combined_df[["D"]][is.na(combined_df[["D"]])] = 0
  return(combined_df)
  
}
bolus_dosis = function(t, dosis, eps=0.1)
{
  return(list(start=t, end=t+eps, rate = dosis / eps))
}

#' Simulates a 1D-SDE on the interval t_start to t_end with given drift and diffusion with noisy observation.
#' @param drift A function which takes time, state and params as input and outputs a scalar.
#' @param diffusion A function which takes time, state and params as input and outputs a scalar.
#' @param params named list with parameters. Should contain at least C0 and sigma_eps
#' @return simulations with shape (n_simulations, t_samples).

simulate_model <- function(drift, diffusion, params, t_start, t_end, n_samples=20, h=0.1, dosis=list( ))
{
  data = get_t_grid(t_start, t_end, n_samples, h, dosis)
  Conc = params$Conc0
  path = c()
  
  for (i in 1:nrow(data))
  { t = data[i, "t"]
    h = data[i + 1, "t"] - data[i, "t"] 
    path = c(path, Conc + rnorm(1, 0, params$sigma_eps**0.5))
    Conc = max(next_step_sde(drift, diffusion,  t, Conc, data[i,"D"], params, h), 0) 
  }
 
  data[["ConcObserved"]] = path
  data[["ConcObserved"]][data[["D_sample"]]] = NaN
  data = data[data[["D_sample"]] | data[["output"]],]
  return(data)
}


  






