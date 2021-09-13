#Dosing and design
get_n_dosis = function(n, d)
{ dosis = list()
  for(i in 0:n)
  {
    dosis[[paste0("element", i)]] = bolus_dosis(i , d, eps=0.01)
  }
  return(dosis)
}

get_samples = function(n, kth_dosis=0, many_samples=FALSE)
{ 
  if(many_samples)
  {
    samples = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35,  0.4, 0.45, 0.5) + kth_dosis
  } else
  {
  samples = c(0.1, 0.2, 0.3, 0.4, 0.5) + kth_dosis
  }

  for(i in 1:n)
{
  samples = c(samples, c(i + 0.1, i - 0.1 ))
}
return(sort(samples))
}
