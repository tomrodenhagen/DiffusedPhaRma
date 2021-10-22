source("diffused_pharma/experiment_utils.R")


#Some helpers for building models

get_sigma_eps_init = function(fixed_sigma_eps)
{
 if(is.null(fixed_sigma_eps))
	{
	  sigma_eps=0.1
	  lower = 0
	  upper = 1
	}
 else
	{
	  sigma_eps = fixed_sigma_eps
	  lower = NULL
	  upper = NULL
	}
 return(list("init"= sigma_eps, "lower" = lower, "upper" = upper))
}
get_H0_model = function(fixed_sigma_eps=NULL, n_retrys=10)
{       sigma_eps_init= get_sigma_eps_init(fixed_sigma_eps)
	model_H0 = build_model(dConc ~ (- CL  *  Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(0)},
                       params=list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_eps"= sigma_eps_init,
				    "sigma_tau"= list("init"= 0.0000001),
				    "Conc0"= list("init"= exp(-5)),
                                   "V"= list("init"= 10, "lower" = 0, "upper" = 100)),
		       n_retrys=n_retrys
			)

	return(model_H0)
}

get_H1_model = function(fixed_sigma_eps=NULL, diffusion_term = "CONSTANT", n_retrys=10)
{       sigma_eps_init= get_sigma_eps_init(fixed_sigma_eps)
	if(diffusion_term == "CONSTANT")
	{
		model_H1 = build_model(dConc ~ (- CL  *  Conc + D / V  ) * dt + sigma_tau * dw1,
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V)},
                       diffusion = function(t, state, u , params){return(0)},
                       params=list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_eps"= sigma_eps_init,
				    "Conc0"= list("init"= exp(-5)),
                                   "V"= list("init"= 10, "lower" = 0, "upper" = 100)),
		        n_retrys=n_retrys
                        )
	}
	if(diffusion_term == "LINEAR")
	{
 	model_H1 = build_model(dConc ~ (- CL - sigma_tau**2 * 0.5  +  D / V / (exp(Conc) )  ) * dt + sigma_tau * dw1, 
                       drift = function(t, state, u, params){return(-params$CL * state + u / params$V  )},
                       diffusion = function(t, state, u , params){return(params$sigma_tau * state  )},
                       params=list("Conc0","CL", "sigma_eps","sigma_tau", "V"),
                       bounds= list("sigma_eps"= sigma_eps_init,
				    "V"= list("init"= 10, "lower" = 0, "upper" = 100),
				    "CL"= list("init"= 1, "lower" = 0.1, "upper" = 10),
				    "Conc0"= list("init"= -5),
				    "sigma_tau"= list("init"= 0.2, "lower" = 0.01, "upper" = 1)
				   ),
		       transf = function(x){return(exp(x))},
                       observation_equation = ConcObserved ~ exp(Conc),
		       n_retrys=n_retrys
                       
	)    
	}
	return(model_H1)
}

 

