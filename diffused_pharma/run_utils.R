library(jsonlite)
library(xtable)



init_experiment_folder = function(config)
{ name = paste(config$name, config$mode, sep="_")
  experiment_folder = file.path(config$base_path, name)
  dir.create(experiment_folder)
  return(experiment_folder)
}
get_run_number = function(run_folder)
{  if(is.null(run_folder))
	{
	 return(0)
         }
   info = read_json(file.path(run_folder, "info.json"))
   return(as.integer(info$run_number))
}

get_last_run_folder = function(experiment_folder)
{
 run_folders = list.dirs(path=experiment_folder, full.names = T, recursive=FALSE)
 if(length(run_folders) == 0)
 {
   return(NULL)
 }
 
 run_numbers = sapply(run_folders, get_run_number)
 return(run_folders[which.max(run_numbers)])
}
is_finished = function(run_folder)
{
 return(list.files(path = run_folder, pattern = "\\.tex$"))
}
init_run_folder= function(experiment_folder, run_number, config)
{
 folder_name = paste("run", run_number, sep="_") 
 folder_path = file.path(experiment_folder, folder_name)
 dir.create(folder_path)
 info_json = list("run_number" = run_number)
 write_json(info_json, 
	  file.path(folder_path, "info.json") )
 write_json(config, file.path(folder_path, "config.json"))
 return(folder_path)
}

prepare_run_folder = function(experiment_folder, config)
{
  last_run_folder = get_last_run_folder(experiment_folder)
  if(config$fresh | is.null(last_run_folder) )
  {
     working_folder = init_run_folder(experiment_folder, get_run_number(last_run_folder) + 1, config)
  }
  else
  {
     working_folder = last_run_folder
  } 
  return(working_folder)
}

get_conf_string = function(res)
{
  return(paste( round(res$rej, 2) , 
	"(", round(res$rej_lower,2), 
	"," , round(res$rej_upper,2),
	")",sep=""))

}
rename_shortcut = function(shortcut)
{
   renaming = list("const_diff" = "K",
		   "linear_diff" = "L",
		   "low_noise" = "Niedriger Noise",
		   "medium_noise" = "Mittlere Noise",
		   "high_noise" = "Hoher Noise",
		   "measure_first_cycle" = "Messung nach 1. Dosis"
		  )
   renamed = renaming[[shortcut]]
   if(is.null(renamed))
   {
	return(shortcut)
   }
   else
   {
    return(renamed) 
   }

}
run_scenarios = function(models, designs, parameter_samplings, config)
{
 experiment_folder = init_experiment_folder(config)
 working_folder = prepare_run_folder(experiment_folder, config)
 res_agg = list()
 for( m in models)
 {
  for(s in parameter_samplings)
  {
    for (d in designs)
    { 
      name = paste(m$shortcut, s$shortcut, d$shortcut, sep="+")
      desc = paste(m$desc, s$desc, d$desc,  sep="::")
      scenario = list(name=name,
                      description=desc,
                      model_H0=m$H0,
                      model_H1=m$H1,
                      design = d$design,
                      sample_params=s$sampling,
                      T_statistic=T_statistic,
                      H0_drift=H0_drift,
                      H1_drift=H1_drift)
      res = run_complete_scenario(scenario, working_folder, config$n_simulations_typ1, config$n_simulations_typ2, config$n_samples, alpha = config$alpha, h=config$h, config$fresh)
      eval_res_typ1 = eval_simulation(res$type1)
      eval_res_typ2 = eval_simulation(res$type2)
      row = list("Modell" = rename_shortcut(m$shortcut), 
		 "Parameter" = rename_shortcut(s$shortcut), 
		 "Design" = rename_shortcut(d$shortcut),
		 "Emp. Typ 1 Fehler" =get_conf_string(eval_res_typ1),
		  "Emp. Power" =  get_conf_string(eval_res_typ2)
		)
      res_agg[[name]] = row
      
    }
  }
 }
 df = do.call(rbind, res_agg)
 rownames(df) = NULL
 write.csv(x=df, file=paste(working_folder, "/res.csv", sep="") )
 print(xtable(df, type = "latex"),
      include.rownames = FALSE,
      floating=FALSE,
       file = file.path(working_folder, paste(config$name, "res.tex", sep="_")) )
}	








