#Visualize scenarios
source("diffused_pharma/experiment_utils.R")
source("diffused_pharma/run_utils.R")
source("scripts/fixtures/models.R")
source("scripts/fixtures/parameters.R")
source("scripts/fixtures/design.R")
source("scripts/config.R")


options <- commandArgs(trailingOnly = TRUE)
config = load_config(options)

#Some script params
test = config$test
n_days = 2


#Some base arguments
low_noise = 0.01
medium_noise = 0.05
high_noise = 0.1

parameter_no_noise = get_parameter_sampling(0)()

#Make CL a bit lower to make it clearer

parameter_no_noise[["CL"]] = 1
drift =function(t, state, u, params){return(-params$CL * state + u / params$V) }
diff =function(t, state, u, params){return(0) }




vis_sigma_eps = function(sigma_eps)
{ dosis = get_n_dosis(n_days, 100)
  design= list(t_start=0, t_end = n_days, n_samples=get_samples(n_days, 0), dosis = dosis ) 
  data_obs = simulate_model(drift, diff, parameter_no_noise, design$t_start, design$t_end , n_samples = design$n_samples, dosis = design$dosis, h=0.01 )
  plot(data_obs[["t"]], data_obs[["ConcObserved"]],
       main= paste("Measure Error for", sigma_eps),
       ylab="Drug concentration in mg /l",
       xlab="time in days", pch=19,col="blue",cex=1.5,type = "n") 
  data_unobs = simulate_model(drift, diff, parameter_no_noise, design$t_start, design$t_end , n_samples = 300, dosis = design$dosis, h=0.01 )
  
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green",lwd = 3
      )
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]] + qnorm(0.95, 0 , sigma_eps**0.5), pch=1, col="red",lwd = 1
  )
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]] - qnorm(0.95, 0 , sigma_eps**0.5), pch=1, col="red",lwd = 1
  )
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green",lwd = 1
  )
}
png(file="./plots/scenario_vis.png", width=5, height=5,
    units="in", res=300, pointsize=8)

par(mfrow=c(3,3), pty="s")

vis_sigma_eps(low_noise)
vis_sigma_eps(medium_noise)
vis_sigma_eps(high_noise)
#dev.off()
#png(file="./plots/diff_cycles_vis.png", width=5, height=5,
#    units="in", res=300, pointsize=8)

#par(mfrow=c(1,3), pty="s")
n_days = 15
vis_different_cycles = function(cycle)
{ dosis = get_n_dosis(n_days, 250)

  design= list(t_start=0, t_end = n_days, n_samples=get_samples(n_days, cycle), dosis = dosis ) 
  data_obs = simulate_model(drift, diff, parameter, design$t_start, design$t_end , n_samples = design$n_samples, dosis = design$dosis, h=0.01 )
  plot(data_obs[["t"]], data_obs[["ConcObserved"]],
       main= paste("Measure in interval", cycle + 1),
       ylab="Drug concentration in mg /l",
       xlab="time in days", pch=19,col="blue",cex=1.5,type = "n") 
  data_unobs = simulate_model(drift, diff, parameter_no_noise, design$t_start, design$t_end , n_samples = 300, dosis = design$dosis, h=0.01 )
   
  lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green",lwd = 3)
  points(data_obs[["t"]], data_obs[["ConcObserved"]], pch=1, col="blue",lwd = 3)
}
vis_different_cycles(0)
vis_different_cycles(3)
vis_different_cycles(6)
#dev.off()

#png(file="./plots/diff_measures_vis.png", width=5, height=5,
#    units="in", res=300, pointsize=8)

#par(mfrow=c(1,3), pty="s")
n_days = 15
vis_different_measures = function(sampling, name)
{ dosis = get_n_dosis(n_days, 250)

design= list(t_start=0, t_end = n_days, n_samples=get_samples(n_days, 3,sampling=sampling), dosis = dosis ) 
data_obs = simulate_model(drift, diff, parameter, design$t_start, design$t_end , n_samples = design$n_samples, dosis = design$dosis, h=0.01 )
plot(data_obs[["t"]], data_obs[["ConcObserved"]],
     main= paste("Measure", name, "number of samples"),
     ylab="Drug concentration in mg /l",
     xlab="time in days", pch=19,col="blue",cex=1.5,type = "n") 
data_unobs = simulate_model(drift, diff, parameter_no_noise, design$t_start, design$t_end , n_samples = 300, dosis = design$dosis, h=0.01 )

lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch=1, col="green",lwd = 3)
points(data_obs[["t"]], data_obs[["ConcObserved"]], pch=1, col="blue",lwd = 3)
}
vis_different_measures("FEW", "few")
vis_different_measures("NORMAL", "normal")
vis_different_measures("MANY", "many")
dev.off()