source("../TrueSkill.R")

# Code 1
mu = 0.0; sigma = 6.0; beta = 1.0; gamma = 0.0; p_draw = 0.0
a1 = Player(Gaussian(mu, sigma), beta, gamma); a2 = Player(); a3 = Player(); a4 = Player()

# Code 2
team_a = c(a1, a2)
team_b = c(a3, a4)
teams = list(team_a, team_b)
g = Game(teams)
# g = Game(teams, c(0,0))

# Code 3
lhs = g$likelihoods
ev = g$evidence
ev = round(ev, 3)
print(ev)

# Code 4
pos = g$posteriors()
print(pos[[1]][1])
print(lhs[[1]][[1]] * p1$prior)

# Code 5
ta = c(a1)
tb = c(a2, a3)
tc = c(a4)
teams = list(ta, tb, tc)
result = c(1, 0, 0)
g = Game(teams, result, p_draw=0.25)

# Code 6
c1 = list(c("a"),c("b"))
c2 = list(c("b"),c("c"))
c3 = list(c("c"),c("a"))
composition = list(c1,c2,c3)
h = History(composition, gamma=0.0)

# Code 7
lc = h$learning_curves()
lc_print(lc$a)
lc_print(lc$b)

# Code 8
h$convergence()
lc = h$learning_curves()
lc_print(lc$a)
lc_print(lc$b)


# Skill evolution

N = 1000
skill <- function(experience, middle, maximum, slope){
    return(maximum/(1+exp(slope*(-experience+middle))))
}
target = skill(seq(N), 500, 2, 0.0075)
opponents = rnorm(N,target,0.5)

composition = list(); results = list(); times = c(); priors = hash()
for(i in seq(N)){composition[[i]] = list(c("a"), c(toString(i)))}
for(i in seq(N)){results[[i]] = if(rnorm(1,target[i])>rnorm(1,opponents[i])){c(1,0)}else{c(0,1)} }
for(i in seq(N)){times = c(times,i)}
for(i in seq(N)){priors[[toString(i)]] = Player(Gaussian(opponents[i],0.2))}
    
h = History(composition, results, times, priors, gamma=0.015)
h$convergence()
lc_a = h$learning_curves()$a; mu = c()
for(tp in lc_a){mu = c(mu,tp[[2]]$mu)}

##############

data = read.csv("input/history.csv", header=T)
get.composition = function(x){
    res = list()
    if (x["double"]=="t"){
        res[[1]] = c(x["w1_name"],x["w2_name"])
        res[[2]] = c(x["l1_name"],x["l2_name"])
    }else{
        res[[1]] = c(x["w1_name"])
        res[[2]] = c(x["l1_name"])
    }
    return(res)
}

