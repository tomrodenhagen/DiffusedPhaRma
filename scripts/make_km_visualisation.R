#Visualize scenarios
source("diffused_pharma/experiment_utils.R")
source("diffused_pharma/run_utils.R")
source("scripts/fixtures/models.R")
source("scripts/fixtures/parameters.R")
source("scripts/fixtures/design.R")
source("scripts/config.R")
library(latex2exp)

options <- commandArgs(trailingOnly = TRUE)
config = load_config(options)

#Some script params
test = config$test
n_days = 15


#Some base arguments
low_noise = 0.01
medium_noise = 0.05
high_noise = 0.1

parameter_no_noise = get_parameter_sampling(0)()

#Make CL a bit lower to make it clearer

dosis = get_n_dosis(n_days, 250)
drift =H1_drift =function(t, state, u, params){return(-params$CL * params$Km / (params$Km + state) * state + u / params$V) }
diff =function(t, state, u, params){return(0) }
sampling="NORMAL"
design= list(t_start=0, t_end = n_days, n_samples=get_samples(n_days, 0,sampling=sampling), dosis = dosis ) 
png(file="./plots/km_vis.png", width=5, height=5,
    units="in", res=300, pointsize=8)

par(mfrow=c(2,2))
for(km in c(3,5,7,9))
{
  parameter_no_noise[["Km"]] = km
data_obs = simulate_model(drift, diff, parameter_no_noise, design$t_start, design$t_end , design$n_samples, dosis = design$dosis, h=0.01 ) 
data_unobs = simulate_model(drift, diff, parameter_no_noise, design$t_start, design$t_end , n_samples = 300, dosis = design$dosis, h=0.01 )

plot(data_obs[["t"]], data_obs[["ConcObserved"]],
      main= TeX(paste("K_m=", km,sep="")),
      ylim=c(0, 15),
      ylab="Drug concentration in mg /l",
      xlab="time in days", pch=19,col="blue",cex=1.5)

lines(data_unobs[["t"]], data_unobs[["ConcObserved"]],
     ylab="Drug concentration in mg /l",
     xlab="time in days", pch=19,col="green",cex=1.5)
abline( h=km,lty = "dashed", col="red"      )

}
dev.off()
