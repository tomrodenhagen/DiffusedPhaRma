library(jsonlite)
library(xtable)
library(latex2exp)


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

get_conf_string = function(point,lower, upper,test=NULL,full=FALSE)
{ if(full)
 { s = paste( round(point, 2) ,
        "(", round(lower,2),
        "," , round(upper,2),
        ")",sep="")
    }else
  {s = paste( round(point, 2), "$\\pm$", round(point - lower,2) ,sep="")
   
  }
  if(is.null(test))
  {
   return(s)
  }
  else
  {
   if(test=="NONE")
   {
    return(s)
   }
   if(test=="BIGGER")
   {
     return(paste("\\colorbox{green!30}{", s, "}", sep=""))
   }
   if(test=="SMALLER")
   {
     return(paste("\\colorbox{red!30}{", s, "}", sep=""))
   }
  }

}
rename_shortcut = function(shortcut, invert=FALSE)
{
   renaming = list("const_diff" = "$\\sigma_\\tau dW$",
		   "linear_diff" = "$C \\cdot \\sigma_\\tau dW$",
		   "const_diff_no_retry" = "$\\sigma_\\tau dW$, kein Multistart",
		   "linear_diff_no_retry" = "$C \\cdot \\sigma_\\tau dW$, kein Multistart",
		   "low_noise" = "$\\sigma_\\epsilon^2 = 0.01$",
		   "medium_noise" = "$\\sigma_\\epsilon^2 = 0.05$",
		   "high_noise" = "$\\sigma_\\epsilon^2 = 0.1$",
		   "const_diff_fixed_low_noise" = "$\\sigma_\\tau dW$,Fix. $\\sigma_\\epsilon^2 = 0.01$",
		   "const_diff_fixed_medium_noise" = "$\\sigma_\\tau dW$,Fix. $\\sigma_\\epsilon^2 = 0.05$",
                   "const_diff_fixed_high_noise" = "$\\sigma_\\tau dW$,Fix. $\\sigma_\\epsilon^2 = 0.1$",
		   "measure_first_cycle" = "$k=1, n=6$",
		   "measure_3th_cycle" = "$k=4, n = 6$",
		   "measure_6th_cylce" = "$k=7, n= 6$",
                   "measure_first_cycle_few_samples" = "$k=1, n=3$",
		   "measure_first_cycle_many_samples" = "$k=1, n=11$"
)  
   if(invert)
   {
   for(i in names(renaming))
   {
    renaming[[renaming[[i]]]] = i
    renaming[[i]] = i
   }
   }
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
init_recent_res_folder = function(config)
{   folder_name =  paste("recent_res", config$mode, sep="_")
    folder_path =  file.path(config$base_path, folder_name )
    dir.create(folder_path)   
    return(folder_path)

}
init_final_res_folder = function(config)
{   folder_name =  paste("final_res", config$mode, sep="_")
    folder_path =  file.path(config$base_path, folder_name )
    dir.create(folder_path)   
    return(folder_path)

}

special_kms = function()
{
 return(c(5,7,9))
}

plot_Km = function(dfs, path, name)
{
 models = lapply(dfs, function(x){x[["Modell"]]} )
 designs = lapply(dfs, function(x){x[["Design"]]} )
 parameters = lapply(dfs, function(x){x[["Parameter"]]} )
 for(m in unique(models))
 { for(p in parameters)
   {df_filtered = dfs[(models==m) & (parameters==p)]
    if(length(df_filtered) > 1)
     {plot_dfs(df_filtered, file.path(path, paste(name, "unfixed",
						  rename_shortcut( m,invert=TRUE),rename_shortcut(p, invert=TRUE) 
						  ,"per_km.png", sep="_")), "Design", m, p)}
   }
 }
 

 for(m in unique(models))
 { for(p in designs)
   {df_filtered = dfs[(models==m) & (designs==p)]
    if(length(df_filtered) > 1)
     {plot_dfs(df_filtered, file.path(path, paste(name, "fixed", rename_shortcut( m,invert=TRUE),rename_shortcut(p, invert=TRUE)
						  ,"per_km.png", sep="_")), "Parameter", m, p )}
   }
 }
 for(m in unique(parameters))
 { for(p in designs)
   {df_filtered = dfs[(parameters==m) & (designs==p)]
    if(length(df_filtered) > 1)
     {plot_dfs(df_filtered, file.path(path, paste(name, "fixed",rename_shortcut( m,invert=TRUE),rename_shortcut(p, invert=TRUE)
						  ,"per_km.png", sep="_")), "Modell", m, p )}
   }
 }

}
plot_dfs = function(dfs,path, unfixed, fixed1,fixed2)
{
 
 png(path)
 plot(NULL, NULL,
       ylab="Test Power",
       xlim=c(2,10),
       ylim=c(0,1.1),
       xlab=TeX("$K_m$"))
 cols = c()
 names = c()
 col_pal = list("red", "green", "blue", "yellow", "brown")
 for(i in seq_along(dfs))
 { df = dfs[[i]]
   col = col_pal[[length(cols)+1]]

   cols = c(cols, col)
   names = c(names,TeX(df[[unfixed]] )) 
   polygon(c(df[["Km"]][,"group"],rev(df[["Km"]][,"group"])), c( df[["Km"]][,"Lower"], rev(df[["Km"]][,"Upper"])),col=col,density=30)
   points(df[["Km"]][,"group"], df[["Km"]][,"Point"], pch=19,col=col,cex=1.5) 
   lines(df[["Km"]][,"group"], df[["Km"]][,"Point"])
   
   }
 
 legend("bottomright", legend=names,pch=15,col=cols)
 dev.off()
}
make_typ2_res_row = function(m,s,d,res,n_sim)
{
  row = list("Modell" = rename_shortcut(m$shortcut), 
	     "Parameter" = rename_shortcut(s$shortcut), 
	     "Design" = rename_shortcut(d$shortcut)
			)
  km_table = res[["per_km"]]
  
  for(i in rownames(km_table))
  { if(km_table[i,"group"] %in% special_kms())
    {
    col_name  = paste("$K_m$=", km_table[i,"group"], sep="")
    row[[col_name]] = list("Point"=km_table[i,"Point"],"Lower"=km_table[i,"Lower"],"Upper"=km_table[i,"Upper"], "N_samples" = km_table[i,"N_samples"])
    }
  row[["Nicht konv. Typ-II"]] = paste(round(res[["not_valid"]] / n_sim* 100,1), "\\%")   
  }
  
  return(row)
				

}
test_for_diff = function(conf1, conf2, alpha=0.05)

{   k1 = conf1$N_samples * conf1$Point
    k2 = conf2$N_samples * conf2$Point
    n1 = conf1$N_samples
    n2 = conf2$N_samples
    less  = fisher.test(rbind(c(k1,n1-k1), c(k2,n2-k2)), alternative="greater") 
    greater  = fisher.test(rbind(c(k1,n1-k1), c(k2,n2-k2)), alternative="less")
    
    if(less$p.value < alpha) 
    {
	return("SMALLER")
    }
    if(greater$p.value < alpha)
    {	   
	 return("BIGGER")
    }
    return("NONE")
   }

add_significance = function(df_typ2, fixed_model="const_diff", fixed_design="measure_first_cycle")
{ 
  fixed_model = rename_shortcut(fixed_model) 
  fixed_design = rename_shortcut(fixed_design)
  df_fixed = df_typ2[df_typ2$Design==fixed_design & df_typ2$Modell == fixed_model,]
  
  if(nrow(df_fixed)==0 | nrow(df_fixed)==nrow(df_typ2))
  { 
    return(make_pretty_res(df_typ2))
  }
  for(col in colnames(df_typ2))
  {
   if(grepl("K_m", col, fixed=TRUE))
     {
	for(row in rownames(df_fixed))
	{
	  km_fixed = df_fixed[row,col][[1]]
	  
	   
	  df_other = df_typ2[df_typ2$Parameter==df_fixed[row,"Parameter"][[1]],]
	  for(row_other in rownames(df_other))
	 { 
	   km_other = df_other[row_other,col][[1]]
	   test_res = test_for_diff(km_fixed, km_other)
	   df_typ2[row_other,col]= get_conf_string(km_other$Point, km_other$Lower, km_other$Upper, test_res)
	   if(df_typ2[row_other,"Design"] ==fixed_design & df_typ2[row_other, "Modell"] == fixed_model)
	   {
	     df_typ2[row_other,col] = paste("\\colorbox{blue!30}{", df_typ2[row_other,col] , "}", sep="")
	   }	   
	 }
	}
     }
  }

return(df_typ2)
}
make_pretty_res= function(df_typ2)
{
 for(col in colnames(df_typ2))
  {
   if(grepl("K_m", col, fixed=TRUE))
   {
    for(row in rownames(df_typ2))
    {
    cell = df_typ2[row,col][[1]]
    df_typ2[row,col] = get_conf_string(as.numeric(cell$Point), cell$Lower, cell$Upper)  
    }
   }
  }
 return(df_typ2)
}

order_df= function(df)
{       for(col in c("Modell", "Design", "Parameter"))
	{
		df[,col] = as.character(df[,col])
	}
	return(df[order(df$Modell,df$Design,df$Parameter),])
}


print_aggregated_results = function(res, recent_res_folder, config)
{
 res_agg_typ1 = res[[1]]
 res_agg_typ2 = res[[2]]
 res_agg_km = res[[3]]
 

 df_typ2 = as.data.frame(do.call(rbind, res_agg_typ2))
 
 df_list = list()
 for(model in unique(df_typ2$Modell))
 {
   df_list[[model]] = add_significance(df_typ2[df_typ2$Modell==model,], fixed_model = model)
 }
 
 df_typ2_design = as.data.frame(do.call(rbind, df_list))
 
 df_list = list()
 for(design in unique(df_typ2$Design))
 {
   df_list[[design]] = add_significance(df_typ2[df_typ2$Design==design,], fixed_design = design)
 }
 df_typ2_model = as.data.frame(do.call(rbind, df_list))
 
 rownames(df_typ2_design) = NULL
 rownames(df_typ2_model) = NULL
 
 df_typ1 = as.data.frame(do.call(rbind, res_agg_typ1))

 
 rownames(df_typ1) = NULL
 
 df_typ2_design = order_df(df_typ2_design)
 df_typ2_model = order_df(df_typ2_model)
 df_typ1 = order_df(df_typ1)
 print(xtable(df_typ2_design, type = "latex"),
      include.rownames = FALSE,
      floating=FALSE,
      file = file.path(recent_res_folder, 
		       paste(config$name, "res_typ2.tex", sep="_")),
       sanitize.colnames.function = identity,
       sanitize.text.function = identity)
 
 print(xtable(df_typ2_model, type = "latex"),
      include.rownames = FALSE,
      floating=FALSE,
      file = file.path(recent_res_folder,
                       paste(config$name, "res_typ2_model.tex", sep="_")),
       sanitize.colnames.function = identity,
       sanitize.text.function = identity)
 print(xtable(df_typ1, type = "latex"),
      include.rownames = FALSE,
      floating=FALSE,
      file = file.path(recent_res_folder, 
		       paste(config$name, "res_typ1.tex", sep="_")),
       sanitize.colnames.function = identity,
       sanitize.text.function = identity  )
 
 plot_Km(res_agg_km, recent_res_folder, config$name) 
		      
}

make_res_path = function(config)
{recent_res_folder = init_recent_res_folder(config)
 return(file.path(recent_res_folder, 
		     paste(config$name, 
			   config$mode, "_res.RDS")) 
	)
}
join_and_compare = function(name1, name2, config1, config2)
{ 
  config1$name = name1
  path1 = make_res_path(config1)
  config2$name = name2
  path2 = make_res_path(config2)
  res1 = readRDS(path1)
  res2 = readRDS(path2)
  res = list()
  res[[1]] = c(res1[[1]], res2[[1]])
  res[[2]] = c(res1[[2]], res2[[2]])
  res[[3]] = c(res1[[3]], res2[[3]]) 
  config1$name = paste("comp", name1, name2, sep="_") 
  print_aggregated_results(res, init_final_res_folder(config1), config1)


}
get_aggregated_results = function(models, designs, parameter_samplings, working_folder, config)
{
 res_agg_typ1 = list()
 res_agg_typ2 = list()
 res_agg_km = list()
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
      res_agg_typ2[[name]] = make_typ2_res_row(m,s,d,eval_res_typ2, config$n_simulations_typ2)
    
      res_agg_typ1[[name]] = list("Modell" = rename_shortcut(m$shortcut), 
		 "Parameter" = rename_shortcut(s$shortcut), 
		 "Design" = rename_shortcut(d$shortcut),
		 "Emp. Typ-I Fehler" =get_conf_string(eval_res_typ1$Point,
						      eval_res_typ1$Lower,
						      eval_res_typ1$Upper),
		  "Nicht konv.Typ-I" = paste(eval_res_typ1[["not_valid"]] /config$n_simulations_typ1 * 100, "\\%")
		 
		)
             
      res_agg_km[[name]] = list("Modell" = rename_shortcut(m$shortcut),
                 "Parameter" = rename_shortcut(s$shortcut),
                 "Design" = rename_shortcut(d$shortcut),
                 "Km" = eval_res_typ2[["per_km"]] ) 
              

            
    }
  }
 }
 return(list(res_agg_typ1, res_agg_typ2, res_agg_km))
}

run_scenarios = function(models, designs, parameter_samplings, config, final_res=FALSE)
{
 experiment_folder = init_experiment_folder(config)
 recent_res_folder = init_recent_res_folder(config)
 working_folder = prepare_run_folder(experiment_folder, config)
 res = get_aggregated_results(models, designs, parameter_samplings, working_folder, config)
 saveRDS(res, make_res_path(config))
 if(final_res)
 {
 print_aggregated_results(res, 
			      init_final_res_folder(config), 
			    config)

 } else
 {
	print_aggregated_results(res, 
			      recent_res_folder, 
			    config)

 }
}	








