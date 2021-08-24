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

build_model = function(equation, drift, diffusion, params, prior=list(), bounds=list(), transf=identity, inv_transf=identity)
{ 
  model <- ctsm$new()
  model$addSystem(equation)

  model$nprior <- length(params)
  model$PriorCorrelationMatrix <- diag(length(params)) 
  dimnames(model$PriorCorrelationMatrix) <- rep(list(params), 2)
  #Conc should be in the equation
  model$addObs( ConcObserved ~ Conc )
  model$addInput(D)
  #sigma_eps should be in the params 
  model$setVariance(ConcObserved ~ sigma_eps)
  set_params(model,params, prior, bounds)
    
  estimate = function(data, n_retrys=5)
    { data[["ConcObserved"]] = inv_transf(data[["ConcObserved"]])
      invisible( capture.output( fit <- model$estimate(data) ) )
      for(k in 1:n_retrys)
      {
        if(fit$info==0)
        {
          break
        }
        if(k==n_retrys)
	{
	 stop("Fitting algorithm didnt converged until last retry.")
	}
        set_params(model,params, prior, bounds)
        invisible( capture.output( fit <- model$estimate(data) ) )

      }
      
      
      return(as.list(fit$xm))
    }
  simulate = function(...)
    {
        data = simulate_model(drift, diffusion, ...)
        data[["ConcObserved"]] = transf(data[["ConcObserved"]])
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
	       "inv_transf" = inv_transf
                  )
  return(myModel)
}
copy_model = function(model)
{
  return(build_model(model$equation, model$drift, 
		     model$diffusion, model$params, 
		     model$prior, model$bounds, 
		     model$transf, model$inv_transf))
}
