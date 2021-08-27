source("diffused_pharma/experiment_utils.R")


model_H1 = build_model(dConc ~ (- CL - sigma_tau**2 * 0.5  +  D / V / exp(Conc)  ) * dt + sigma_tau * dw1, 
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V  )},
                       diffusion = function(t, state, u , params){return(params$sigma_tau * state  )},
                       list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("V"= list("init"= 10, "lower" = 0, "upper" = 100),
                       "Conc0"= list("init"= 0, "lower" = -10, "upper" = 10)),
                       
		       transf = function(x){return(exp(x))},
                       observation_equation = ConcObserved ~ exp(Conc)
                       
)
 


sample_params = function()
{ V = 1
Km = 4 # mg per liter
Vmax = 5
CL = Vmax / Km # mg per day
Conc0 = 0
return(list("sigma_eps"=0.01, "sigma_tau"=1, "CL"=CL, "Conc0"=Conc0, "Km"= Km, "V"=V))
}
params = sample_params()

dosis = list(bolus_dosis(0.06, 5),bolus_dosis(2.06, 5))
data_H1 =model_H1$simulate(params, 0, 3, 100, h=0.01, dosis)
plot(data_H1[["ConcObserved"]])
est=model_H1$estimate(data_H1)
print(est)
