#Dosing and design
get_n_dosis = function(n, d)
{ dosis = list()
  for(i in 0:n)
  {
    dosis[[paste0("element", i)]] = bolus_dosis(i , d, eps=0.01)
  }
  return(dosis)
}

get_samples = function(n, kth_dosis=0, sampling="NORMAL")
{ 
  if(sampling=="MANY")
  {
    samples = c( 0.15, 0.2, 0.25, 0.3, 0.35,  0.4, 0.45, 0.5,0.55, 0.6, 0.65) + kth_dosis
  } 
  if(sampling=="NORMAL")
  {
  samples = c(0.15, 0.25, 0.35, 0.45, 0.55, 0.65) + kth_dosis
  }
  if(sampling=="FEW")
  {
  samples = c( 0.2, 0.4, 0.6) + kth_dosis
  }
  for(i in 1:n)
{
  samples = c(samples, c(i + 0.1, i - 0.1 ))
}
  samples = c(samples, 0.1)
return(sort(samples))
}
