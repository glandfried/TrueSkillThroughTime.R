source("../TrueSkill.R")
library(microbenchmark)
library(profvis)


print("Code 1")
mu = 0.0; sigma = 6.0; beta = 1.0; gamma = 0.03; p_draw = 0.0

print("Code 2")
a1 = Player(Gaussian(mu, sigma), beta, gamma); a2 = Player(); a3 = Player(); a4 = Player()

print("Code 3")
team_a = c(a1, a2)
team_b = c(a3, a4)
teams = list(team_a, team_b)
g = Game(teams)
# g = Game(teams, c(0,0))
microbenchmark(posteriors(Game(teams)), times=10, unit="s")
#profvis(microbenchmark(posteriors(Game(teams)), times=100, unit="s"))

print("Code 4")
lhs = g@likelihoods
ev = g@evidence
ev = round(ev, 3)

print("Code 5")
pos = posteriors(g)
print(pos[[1]][1])
print(lhs[[1]][[1]] * a1@prior)

print("Code 6")
ta = c(a1)
tb = c(a2, a3)
tc = c(a4)
teams = list(ta, tb, tc)
result = c(1, 0, 0)
g = Game(teams, result, p_draw=0.25)
microbenchmark(posteriors(Game(teams, result, p_draw=0.25)), times=10, unit="s")

#profvis(Game(teams))
#profvis(microbenchmark(posteriors(Game(teams, result, p_draw=0.25)), times=20, unit="s"))

print("Code 7")
c1 = list(c("a"),c("b"))
c2 = list(c("b"),c("c"))
c3 = list(c("c"),c("a"))
composition = list(c1,c2,c3)
h = History(composition, gamma=0.0)
microbenchmark(History(composition, gamma=0.0), times=20, unit="s")

print("Code 8")
lc = h$learning_curves()
lc_print(lc$a)
lc_print(lc$b)

print("Code 9")
h$convergence()
microbenchmark(h$convergence(iterations=1,verbose=F), times=10, unit="s")
lc = h$learning_curves()
lc_print(lc$a)
lc_print(lc$b)

# Skill evolution
print("Code 10")
N = 1000
skill <- function(experience, middle, maximum, slope){
    return(maximum/(1+exp(slope*(-experience+middle))))
}
target = skill(seq(N), 500, 2, 0.0075)
opponents = rnorm(N,target,0.5)

print("Code 11")
composition = list(); results = list(); times = c(); priors = hash()
for(i in seq(N)){composition[[i]] = list(c("a"), c(toString(i)))}
for(i in seq(N)){results[[i]] = if(rnorm(1,target[i])>rnorm(1,opponents[i])){c(1,0)}else{c(0,1)} }
for(i in seq(N)){times = c(times,i)}
for(i in seq(N)){priors[[toString(i)]] = Player(Gaussian(opponents[i],0.2))}

h = History(composition, results, times, priors, gamma=0.015)
#microbenchmark(History(composition, results, times, priors, gamma=0.015), times=1, unit="s")
#microbenchmark(h$convergence(iterations=1), times=1, unit="s")
lc_a = h$learning_curves()$a; mu = c()
for(tp in lc_a){mu = c(mu,tp[[2]]@mu)}
