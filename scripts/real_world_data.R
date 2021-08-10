source("diffused_pharma/experiment_utils.R")
library("stringr")
model_H0 = build_model(dConc ~ (- CL  * Conc + D  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state)},
                       diffusion = function(t, state, u , params){return(0)},
                       list("Conc0","CL", "sigma_eps","sigma_tau"),
                       bounds= list("sigma_tau"= list("init"= 0.0000001))
)

model_H1 = build_model(dConc ~ (- CL  * Conc + D  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state)},
                       diffusion = function(t, state, u , params){return(params$sigma_tau)},
                       list("Conc0","CL", "sigma_eps","sigma_tau"))

T_statistic = function(estimate){return(estimate$sigma_tau)[0]}

design = list(t_start=0, t_end = 10, n_samples=20, dosis = list() )

string_to_list = function(string)
{ 
  string = str_replace_all(string, "\\[", "" )
  string = str_replace_all(string, "\\]", "" )
  string = str_replace_all(string, " ", "" )
  
  return(as.numeric(unlist(strsplit(string, ","))))
}

process_row = function(row)
{ k=6
  time = string_to_list(row[k:k+1,"time"]) 
  print(row[k:k+1,"time"])
  print(row)
  Conc = string_to_list(row[k:k+1,"mean"]) 
 
  plot(time, Conc)
  return(na.omit(data.frame("t"=time, "ConcObserved"=Conc, "D"=rep(0,length(time)))))
  
}
data = read.csv("C:/Users/roden/Downloads/pkdb_data3d0f86d3-6104-4fa9-991e-060ee4c7313a/timecourses.csv")

data = data[data[["substance"]] == "warfarin",]
df = process_row(data)
print(df)
design  = list(t_start = 0, t_end = 400, dosis=list(), n_samples = df[["t"]])
T_statistic = function(estimate){return(estimate$sigma_tau)[0]}
print(run_test(df, model_H0, model_H1, T_statistic=T_statistic, n_simulations=100 ,design=design))

