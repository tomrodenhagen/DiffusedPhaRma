



load_config = function(cmd_args)
{
 if(length(cmd_args)==0)
 {
   s = "TEST"
 }
 else
 { s = cmd_args[1] } 
 local = .Platform$OS.type != "unix"
 if(!local)
 {
   base_path = "/localhome/tr/runs"
 }
 else
 {
   base_path = "C:/Users/roden/Documents/data"
 }
 base_config = list(base_path =  base_path, 
		    fresh = !("CONT" %in% cmd_args),
		    n_days= 15,
		    h = 0.02,
		    alpha = 0.05)
 
 if(s=="TEST")
 {
   base_config$n_simulations_typ1 =  50
   base_config$n_simulations_typ2 = 50
   base_config$n_samples = 50
   base_config$n_retries = 0
 
 }
 if(s=="FULL")
 {
   base_config$n_simulations_typ1 =  100
   base_config$n_simulations_typ2 = 200
   base_config$n_samples = 200
   base_config$n_retries = 4

 }
  if(s=="FULL_EXT")
 {
   base_config$n_simulations_typ1 =  200
   base_config$n_simulations_typ2 = 500
   base_config$n_samples = 500
   base_config$n_retries = 4

 }
 if(s=="DEBUG")
 {
   
   base_config$n_simulations_typ1 =  25
   base_config$n_simulations_typ2 = 25
   base_config$n_samples = 25
   base_config$n_retries = 0
 
 }
 base_config$mode = s
 
 return(base_config)
}
