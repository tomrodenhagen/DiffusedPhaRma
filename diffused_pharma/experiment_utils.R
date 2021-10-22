source("diffused_pharma/fit.R")
source("diffused_pharma/sim.R")
library(parallel)
library(mctools)
library(jsonlite)
DEBUG = FALSE
linux = .Platform$OS.type == "unix"
if (linux) {
    library(unixtools)
}
run_test = function(data, model_H0, model_H1, T_statistic, n_simulations, design,
    alpha = 0.05, h = 0.1) {
    estimated_params_H0 = model_H0$estimate(data)
    H1_fit = model_H1$estimate(data, full_info = TRUE)
    if (is.null(estimated_params_H0) | is.null(H1_fit)) {
        return(NULL)
    }
    estimated_params_H1 = H1_fit$res
    T_samples = c()
    generate_T_sample = function(dummy) {
        data_sim = model_H0$simulate(estimated_params_H0, design$t_start, design$t_end,
            design$n_samples, h, design$dosis)
        estimate_sim = model_H1$estimate(data_sim)
        if (!is.null(estimate_sim)) {
            return(T_statistic(estimate_sim))
        } else {
            return(NA)
        }
    }
    T_samples <- mclapply(1:n_simulations, generate_T_sample, mc.cores = 1)
    T_samples = unlist(T_samples, use.names = FALSE)
    emp_quantile = unname(quantile(na.omit(T_samples), 1 - alpha))
    return(list(rejected = emp_quantile < T_statistic(estimated_params_H1), H1_lh = H1_fit$f,
        all_samples = na.omit(T_samples), t_statistic = T_statistic(estimated_params_H1),
        quantile = emp_quantile, na_samples = sum(is.na(T_samples))))
}
visualize_setting = function(drift, diffusion, model_H0, model_H1, T_statistic, sample_params,
    design, h, path, draw_Km = FALSE) {
    sampled_params = sample_params()
    data_obs = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end,
        design$n_samples, h = h, design$dosis)
    sampled_params[["sigma_eps"]] = 0
    data_unobs = simulate_model(drift, diffusion, sampled_params, design$t_start,
        design$t_end, 1000, h = h, design$dosis)
    estimated_params_H0 = model_H0$estimate(data_obs)
    estimated_params_H1 = model_H1$estimate(data_obs)
    if (is.null(estimated_params_H0) | is.null(estimated_params_H1)) {
        return(NULL)
    }
    path_only_design = file.path(path, "only_design.png")
    jpeg(file = path_only_design)
    plot(data_obs[["t"]], data_obs[["ConcObserved"]], main = "Sample experiment",
        ylab = "Drug concentration in mg /l", xlab = "time in days", pch = 19, col = "blue",
        cex = 1.5)
    lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch = 1, col = "green")
    legend(x = "bottomright", legend = c("Observed points", "True"), col = c("blue",
        "green"), lty = 1:2, cex = 1, pch = 19)
    dev.off()
    pathH0 = file.path(path, "visH0.png")
    jpeg(file = pathH0)
    plot(data_obs[["t"]], data_obs[["ConcObserved"]], main = "Sample experiment",
        ylab = "Drug concentration in mg /l", xlab = "time in days", pch = 19, col = "blue",
        cex = 1.5)
    lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch = 1, col = "green")
    data_H0 = model_H0$simulate(estimated_params_H0, design$t_start, design$t_end,
        300, h = h, design$dosis)
    lines(data_H0[["t"]], data_H0[["ConcObserved"]], pch = 2, col = "red")
    legend(x = "bottomright", legend = c("Observed points", "True", "Predicted from Model H0",
        "Km"), col = c("blue", "green", "red", "brown"), lty = 1:2, cex = 1, pch = 19)
    if (draw_Km) {
        abline(h = sampled_params$Km, col = "brown")
    }
    dev.off()
    pathH1 = file.path(path, "visH1.png")
    jpeg(file = pathH1)
    plot(data_obs[["t"]], data_obs[["ConcObserved"]], main = "Sample experiment",
        ylab = "Drug concentration in mg /l", xlab = "time in days", pch = 19, col = "blue",
        cex = 1.5)
    lines(data_unobs[["t"]], data_unobs[["ConcObserved"]], pch = 1, col = "green")
    data_H1 = model_H1$simulate(estimated_params_H1, design$t_start, design$t_end,
        300, h = h, design$dosis)
    lines(data_H1[["t"]], data_H1[["ConcObserved"]], pch = 3, col = "blue")
    if (draw_Km) {
        abline(h = sampled_params$Km, col = "brown")
    }
    legend(x = "bottomright", legend = c("Observed points", "True", "Predicted from Model H1",
        "Km"), col = c("blue", "green", "blue", "brown"), lty = 1:2, cex = 1)
    dev.off()
}
library(progress)
run_simulation_study = function(drift, diffusion, model_H0, model_H1, T_statistic,
    sample_params, design, n_param_samples, n_simulations, alpha, log_path, h = 0.1) {
    rec <- as.data.frame(matrix(0, ncol = length(sample_params()), nrow = n_param_samples))
    names(rec) = names(sample_params())
    res_vec = c()
    pb = progress_bar$new(total = n_param_samples)
    linux = .Platform$OS.type == "unix"
    if (linux) {
        numCores = detectCores()
    } else {
        numCores <- 1
    }
    run_test_wrapper = function(dummy) {
        if (linux) {
            tmp_path = file.path("/tmp", paste("tmp", dummy, sep = "_"))
            dir.create(tmp_path, showWarnings = FALSE)
            set.tempdir(tmp_path)
        }
        sampled_params = sample_params()
        data = simulate_model(drift, diffusion, sampled_params, design$t_start, design$t_end,
            design$n_samples, h = h, design$dosis)
        # Data should look like real world data
        model_H0_copy = copy_model(model_H0)
        model_H1_copy = copy_model(model_H1)
        res = run_test(data, model_H0_copy, model_H1_copy, T_statistic, n_simulations,
            design, alpha, h)
        if (linux) {
            unlink(tmp_path, recursive = TRUE)
        }
        return(list(res = res, params = unlist(sampled_params, use.names = FALSE)))
    }
    res_list = mcMap(1:n_param_samples, run_test_wrapper, mc.cores = numCores)
    res_vec = rep(NA, length(res_list))
    for (i in seq_along(res_list)) {
        rec[i, ] = res_list[[i]]$params
        if (!is.null(res_list[[i]]$res)) {
            res_vec[[i]] = res_list[[i]]$res$rejected
        }
    }
    rec[["test_rejected"]] = as.integer(res_vec)
    return(rec)
}
library(matrixStats)
binconf = function(k, n, alpha = 0.05) {
    # Use normal distributin approximation
    p = k/n
    b = (p * (1 - p)/n)^0.5 * qnorm(1 - alpha/2)
    return(list(Point = p, Upper = p + b, Lower = p - b, N_samples = n))
}
conf_per_group = function(x, g) {
    ks = aggregate(x, by = list(g), FUN = sum)
    ns = aggregate(x, by = list(g), FUN = length)
    res = list()
    for (i in rownames(ns)) {
        conf = binconf(ks[i, "x"], ns[i, "x"])
        conf[["group"]] = ks[i, 1]
        res[[i]] = conf
    }
    res = as.data.frame(do.call(rbind, res))
    res$group = as.numeric(res$group)
    res$Point = as.numeric(res$Point)
    res$Lower = as.numeric(res$Lower)
    res$Upper = as.numeric(res$Upper)
    res$N_samples = as.numeric(res$N_samples)
    return(res)
}
eval_simulation = function(rec_with_na, path = NULL, name = NULL) {
    n_not_valid = sum(is.na(rec_with_na$test_rejected))
    rec = rec_with_na[!is.na(rec_with_na$test_rejected), ]
    conf = binconf(sum(rec$test_rejected), length(rec$test_rejected))
    res_evaluated = list(Point = mean(rec$test_rejected), Upper = conf$Upper, Lower = conf$Lower,
        not_valid = n_not_valid)
    if (!is.null(path)) {
        write(toJSON(res_evaluated), file.path(path, paste(name, "eval.csv", sep = "_")))
    }
    conf_per_Km = conf_per_group(rec$test_rejected, rec$Km)
    res_evaluated[["per_km"]] = conf_per_Km
    return(res_evaluated)
}
set_up_scenario_folder = function(scenario, path, fresh) {
    scenario_folder = file.path(path, scenario$name)
    typ1_folder = file.path(scenario_folder, "typ1")
    typ2_folder = file.path(scenario_folder, "typ2")
    if (dir.exists(scenario_folder) & fresh) {
        unlink(scenario_folder, rec = TRUE)
    }
    if (!dir.exists(scenario_folder)) {
        dir.create(scenario_folder)
        dir.create(typ1_folder)
        dir.create(typ2_folder)
    }
    return(list(scenario_folder = scenario_folder, typ1_folder = typ1_folder, typ2_folder = typ2_folder))
}
run_complete_scenario = function(scenario, path = "C:/Users/roden/Dropbox/Masterarbeit",
    n_tests_typ1 = 500, n_tests_typ2 = 500, n_samples = 500, alpha = 0.05, h = 0.02,
    fresh = TRUE) {
    folders = set_up_scenario_folder(scenario, path, fresh)
    cat(scenario$description, file = file.path(folders$scenario_folder, "description.txt"),
        sep = "\n")
    # some standart parameters, we might change that at a later point
    diffusion = function(t, state, u, params) {
        return(0)
    }
    # Visualize
    if (fresh | !file.exists(file.path(folders$typ1_folder, "visH0.png"))) {
        visualize_setting(scenario$H0_drift, diffusion, scenario$model_H0, scenario$model_H1,
            scenario$T_statistic, scenario$sample_params, scenario$design, h, folders$typ1_folder,
            )
    }
    if (fresh | !file.exists(file.path(folders$typ1_folder, "visH1.png"))) {
        visualize_setting(scenario$H1_drift, diffusion, scenario$model_H0, scenario$model_H1,
            scenario$T_statistic, scenario$sample_params, scenario$design, h, folders$typ2_folder,
            draw_Km = TRUE)
    }
    # Type 1 Error
    res_path = file.path(folders$typ1_folder, "complete_res.csv")
    if (file.exists(res_path) & !fresh) {
        res_type_1 = readRDS(file = res_path)
    } else {
        print("Run type1 testing...")
        res_type_1 = run_simulation_study(scenario$H0_drift, diffusion, scenario$model_H0,
            scenario$model_H1, scenario$T_statistic, scenario$sample_params, scenario$design,
            n_tests_typ1, n_samples, alpha, folders$typ2_folder, h = h)
        saveRDS(res_type_1, file = res_path)
    }
    eval_simulation(res_type_1, folders$scenario_folder, "typ1_res")
    # Type 2 error
    res_path = file.path(folders$typ2_folder, "complete_res.csv")
    if (file.exists(res_path) & !fresh) {
        res_type_2 = readRDS(file = res_path)
    } else {
        print("Running type 2 testing...")
        res_type_2 = run_simulation_study(scenario$H1_drift, diffusion, scenario$model_H0,
            scenario$model_H1, scenario$T_statistic, scenario$sample_params, scenario$design,
            n_tests_typ2, n_samples, alpha, folders$typ2_folder, h = h)
        saveRDS(res_type_2, file = res_path)
    }
    eval_simulation(res_type_2, folders$scenario_folder, "typ2_res")
    return(list(type1 = res_type_1, type2 = res_type_2))
}
