# Visualize scenarios
source("diffused_pharma/experiment_utils.R")
source("diffused_pharma/run_utils.R")
source("scripts/fixtures/models.R")
source("scripts/fixtures/parameters.R")
source("scripts/fixtures/design.R")
source("scripts/config.R")
drift = function(t, state) {
    return(-state)
}
diff = function(t, state) {
    return(0.1)
}
simulate = function(x_start, drift, diff, h, t_end) {
    t = 0
    x = x_start
    ts = c(t)
    xs = c(x_start)
    while (t < t_end) {
        x = x + h * drift(t, x) + h^0.5 * rnorm(1, 0, 1) * diff(t, x)
        t = t + h
        xs = c(xs, x)
        ts = c(ts, t)
    }
    return(list(t = ts, x = xs))
}
png(file = "./plots/vis_simulation.png", width = 9, height = 3, units = "in", res = 300,
    pointsize = 8)
par(mfrow = c(1, 3), pty = "s")
data_h05 = simulate(2, drift, diff, 0.5, 10)
data_h005 = simulate(2, drift, diff, 0.05, 10)
data_h0005 = simulate(2, drift, diff, 0.005, 10)
plot(data_h01$t, data_h01$x, type = "l", xlab = "t", ylab = "X", main = "h=0.5",
    cex.main = 2.5)
plot(data_h001$t, data_h001$x, type = "l", xlab = "t", ylab = "X", main = "h=0.05",
    cex.main = 2.5)
plot(data_h0001$t, data_h0001$x, type = "l", xlab = "t", ylab = "X", main = "h=0.005",
    cex.main = 2.5)
dev.off()
