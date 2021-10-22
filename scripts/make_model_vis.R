source("diffused_pharma/sim.R")
source("diffused_pharma/experiment_utils.R")
H0_drift = function(t, state, u, params) {
    return(-params$CL * state + u/params$V)
}
H1_drift = function(t, state, u, params) {
    return(-params$CL * params$Km/(params$Km + state) * state + u/params$V)
}
sample_params = function() {
    V = 50
    Km = 4  # mg per liter
    Vmax = 5
    CL = Vmax/Km  # mg per day
    Conc0 = 0
    return(list(sigma_eps = 0, sigma_tau = 0, CL = CL, Conc0 = Conc0, Km = Km, V = V))
}
params = sample_params()
dosis = list(bolus_dosis(0, 500, eps = 0.01))
data_H0 = simulate_model(H0_drift, function(...) {
    0
}, params, 0, 5, 500, h = 0.01, dosis)
data_H1 = simulate_model(H1_drift, function(...) {
    0
}, params, 0, 5, 500, h = 0.01, dosis)
print(data_obs)
plot(data_H1[["t"]], data_H1[["ConcObserved"]], type = "l", col = "green", xlab = "Time",
    ylab = "Concentration")
lines(data_H0[["t"]], data_H0[["ConcObserved"]], col = "red")
abline(h = 4, lty = "dashed", col = "blue")
legend(x = "topright", legend = c("Linear Model", "Nonlinear Model", "Km"), col = c("red",
    "green", "blue"), lty = 1:2, cex = 1, pch = 19)
plot(data_H1[["t"]], log(data_H1[["ConcObserved"]]), type = "l", col = "green", xlab = "Time",
    ylab = "Log of Concentration")
lines(data_H0[["t"]], log(data_H0[["ConcObserved"]]), col = "red")
legend(x = "topright", legend = c("Linear Model", "Nonlinear Model", "Km"), col = c("red",
    "green", "blue"), lty = 1:2, cex = 1, pch = 19)
abline(h = log(4), lty = "dashed", col = "blue")
# Choice of parameters
