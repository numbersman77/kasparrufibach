# script to illustrate non-proportional hazards in biomarkers subgroups

# biomarker negatives
m1 <- 4
hr1 <- 0.8
m2 <- m1 / hr1
lam1 <- log(2) / m1
lam2 <- log(2) / m2

# biomarker positives
m3 <- 8
hr2 <- 0.5
m4 <- m3 / hr2
lam3 <- log(2) / m3
lam4 <- log(2) / m4

xs <- seq(0, 100, by = 0.01)
S1 <- 1 - pexp(xs, rate = lam1)
S2 <- 1 - pexp(xs, rate = lam2)
S3 <- 1 - pexp(xs, rate = lam3)
S4 <- 1 - pexp(xs, rate = lam4)

# plot survival functions by biomarkers subgroup and treatment
par(mfrow = c(1, 3), las = 1, cex = 1)
plot(0, 0, type = "n", xlim = c(0, 60), ylim = c(0, 1), xlab = "time", ylab = "event-free probability", main = "survival functions in biomarkers subgroups")

lines(xs, S1, type = "l", col = 2, lty = 3, lwd = 3)
lines(xs, S2, type = "l", col = 3, lty = 3, lwd = 3)

lines(xs, S3, type = "l", col = 2, lty = 1, lwd = 3)
lines(xs, S4, type = "l", col = 3, lty = 1, lwd = 3)

legend("topleft", c("bio-, control", "bio-, treatment", "bio+, control", "bio+, treatment"), bty = "n", 
       col = c(2, 3, 2, 3), lty = c(3, 3, 1, 1), seg.len = 1, lwd = 3)

# marginal survival functions
p <- 1 / 2
S5 <- p * S1 + (1 - p) * S3
S6 <- p * S2 + (1 - p) * S4

plot(0, 0, type = "n", xlim = c(0, 60), ylim = c(0, 1), xlab = "time", ylab = "event-free probability", main = "marginal survival functions by treatment group")
lines(xs, S5, type = "l", col = 2, lty = 1, lwd = 3)
lines(xs, S6, type = "l", col = 3, lty = 1, lwd = 3)

# compute hazards using h(t) = dS(t) / dt * (S(t)) ^ (-1)
h5 <- (p * lam1 * S1 + (1 - p) * lam3 * S3) / (p * S1 + (1 - p) * S3)
h6 <- (p * lam2 * S2 + (1 - p) * lam4 * S4) / (p * S2 + (1 - p) * S4)

# plot hazard ratio
plot(xs, h6 / h5, type = "l", xlab = "time", ylab = "hazard ratio", lwd = 3, col = 2, main = "hazard ratio", xlim = c(0, 60))



