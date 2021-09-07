#Some helper functions
source("diffused_pharma/sim.R")
library(ctsmr)

set_params = function(model, params, prior, bounds)
{
  for (name in params)
  { 
    if(!is.null(bounds[[name]]))
    { 
      init = bounds[[name]][["init"]]
      lower = bounds[[name]][["lower"]]
      upper = bounds[[name]][["upper"]]
    }
    else
    { 
      
      lower = 0
      upper = 10
    }
    if(!is.null(prior[[name]]))
    {
      mu = prior[[name]][["mean"]]
      sigma = prior[[name]][["sd"]]
      args=list( c(init=mu, lower, upper, psd=sigma ) )
      names(args) = list(name)
      do.call(model$setParameter, args)
      
    }
    else
      
    { if(is.null(lower))
    {
      args =list( c(init=init) )
    }
      else
      {
        args =list( c(init = runif(1, lower , upper), lower, upper ) )
      }
      names(args) = list(name)
      do.call(model$setParameter, args)
      
    }
  }
}

build_model = function(equation, drift, diffusion, params, prior=list(), bounds=list(), observation_equation=ConcObserved ~ Conc, transf=identity)
{ 
  model <- ctsm$new()
  model$addSystem(equation)

  model$nprior <- length(params)
  model$PriorCorrelationMatrix <- diag(length(params)) 
  dimnames(model$PriorCorrelationMatrix) <- rep(list(params), 2)
  #Conc should be in the equation
  model$addObs(observation_equation )
  model$addInput(D)
  #sigma_eps should be in the params 
  model$setVariance(ConcObserved ~ sigma_eps)
  set_params(model,params, prior, bounds)
    
  estimate = function(data, n_retrys=10, full_info=FALSE)
    { 
      
      best_fit = NULL
      for(k in 1:n_retrys)
      { invisible( capture.output( fit <- model$estimate(data) ) ) 
        if(is.null(best_fit))
	{
	 if(fit$info==0)
	  {
	   best_fit=fit
	  }
	}
        else
        {
          if(fit$info==0)
	{
	 
	 if(best_fit$f > fit$f)
	  {
	  best_fit=fit
	  }
		 
        }
     
       
        set_params(model,params, prior, bounds)
       }
      
      }
      if(is.null(best_fit)) 
	{
	 warning("Fitting algorithm didnt converged until last retry.")
	 return(NULL)
	}
      res = as.list(best_fit$xm)
      res[["Conc0"]] = transf(res[["Conc0"]])
      if(full_info)
      {       best_fit[["res"]] = res
	      return(best_fit)}
      else
      {
      return(res)
      }
    }
  simulate = function(...)
    {
        data = simulate_model(drift, diffusion, ...)
	return(data)
      
    } 

 
  
  myModel = list("cstmModel" = model,
               "drift"=drift,
               "diffusion" = diffusion,
               "estimate"  = estimate,
               "simulate"  = simulate,
	       "equation" = equation,
	       "params" = params,
	       "prior" = prior,
	       "bounds" = bounds,
	       "transf"= transf,
	       "observation_equation" = observation_equation
                  )
  return(myModel)
}
copy_model = function(model)
{
  return(build_model(model$equation, model$drift, 
		     model$diffusion, model$params, 
		     model$prior, model$bounds, 
		      model$observation_equation, model$transf))
}
