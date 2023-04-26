library(TrueSkillThroughTime)

# First solve your own example. Here is a dummy one.
agents <- c("a", "b", "c", "d", "e")
composition <- list()
for (i in 1:500) {
  who = sample(agents, 2)
  composition[[i]] <- list(list(who[1]), list(who[2]))
}
h <- History(composition = composition, gamma = 0.03, sigma = 1.0)
h$convergence(iterations=6)

# Then plot some learning curves
lc <- h$learning_curves()
colors <- c(rgb(0.2,0.2,0.8), rgb(0.2,0.8,0.2), rgb(0.8,0.2,0.2))
colors_alpha <- c(rgb(0.2,0.2,0.8,0.2), rgb(0.2,0.8,0.2,0.2), rgb(0.8,0.2,0.2,0.2))
plot(0,0, xlim = c(0, 500), ylim = c(-1, 1), xlab = "t", ylab = "skill", type = "n")
for (i in 1:3) {
  agent <- agents[i]
  t <- c(); mu <- c(); sigma <- c()
  for(x in lc[[agent]]){
    t <- c(t, x$t )
    mu <- c(mu, x$N@mu)
    sigma <- c(sigma, x$N@sigma)
  }
  lines(t, mu, col = colors[i], lwd = 2, type = "l")
  polygon(c(t, rev(t)), c(mu + sigma, rev(mu - sigma)), col = colors_alpha[i], border = NA)
}
legend("topright", legend = agents[1:3], col = colors, lwd = 2)
