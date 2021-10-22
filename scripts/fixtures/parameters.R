# Choice of parameters
get_parameter_sampling = function(sigma_eps, fixed_Km = NULL) {
    sample_params = function() {
        V = 50
        if (is.null(fixed_Km)) {
            Km = sample(c(3, 4, 5, 6, 7, 8, 9), 1)
        } else {
            Km = fixed_Km
        }
        Vmax = 7
        CL = Vmax/4  # mg per day
        Conc0 = 0
        return(list(sigma_eps = sigma_eps, sigma_tau = 0, CL = CL, Conc0 = Conc0,
            Km = Km, V = V))
    }
    return(sample_params)
}
