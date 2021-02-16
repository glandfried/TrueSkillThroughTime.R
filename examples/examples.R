source("../TrueSkill.R")

# Code 1
mu = 0.0; sigma = 6.0; beta = 1.0; gamma = 0.0
p1 = Player(Gaussian(mu, sigma), beta, gamma); p2 = Player(Gaussian(mu, sigma), beta, gamma);
p3 = Player(Gaussian(mu, sigma), beta, gamma); p4 = Player(Gaussian(mu, sigma), beta, gamma);

# Code 2
team_a = c(p1, p2)
team_b = c(p3, p4)
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
ta = c(p1)
tb = c(p2, p3)
tc = c(p4)
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
lc$a
lc$b
