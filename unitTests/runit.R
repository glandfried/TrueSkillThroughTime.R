source("TrueSkill.R")
if (!require("RUnit", quietly = TRUE)) {
  stop("Package Runit is not found.") 
}

test_gaussian_init = function() {
    checkEquals(Gaussian(0,1)$isapprox(N01), T)
    checkEquals(Gaussian(0,Inf)$mu, 0)
    checkEquals(Gaussian(0,Inf)$sigma == Inf, T)
    checkEquals(Gaussian(0,0)$isapprox(N00), T)
}
test_ppf = function(){
    checkEquals(ppf(0.3,0, 1),-0.52440044)
    checkEquals(ppf(0.3,2,3),0.42679866)
}
test_cdf = function(){
    checkEquals(cdf(0.3,0,1),0.617911409)
    checkEquals(cdf(0.3,2,3),0.28547031)
}
test_pdf = function(){
    checkEquals(pdf(0.3,0,1),0.38138781)
    checkEquals(pdf(-0.3,0,1),0.38138781)
    checkEquals(pdf(0.3,2,3),0.11325579)
}
test_compute_margin = function(){
    checkEquals(compute_margin(0.25,sqrt(2)*25.0/6),1.8776005988)
    checkEquals(compute_margin(0.25,sqrt(3)*25.0/6),2.29958170)
    checkEquals(compute_margin(0.0,sqrt(3)*25.0/6),2.7134875810e-07)
}
test_trunc = function(){
    checkEquals(trunc(0,1,0.,FALSE),Gaussian(0.797884536,0.6028103066) )
    checkEquals(trunc(0.,sqrt(2)*(25/6),1.8776005988,TRUE), Gaussian(0.0,1.0767055), tolerance = 1e-6)
    checkEquals(trunc(12,sqrt(2)*(25/6),1.8776005988,TRUE),Gaussian( 0.3900995,1.0343979),tolerance = 1e-6)
}
test_gaussian = function(){
    N = Gaussian(25.0, 25.0/3); M = Gaussian(0.0, 1.0)
    checkTrue((M/N)$isapprox(Gaussian(-0.365, 1.007),1e-3))
    checkTrue((M*N)$isapprox(Gaussian(0.355, 0.993),1e-3))
    checkTrue((M+N)$isapprox(Gaussian(25.000, 8.393),1e-3))
    checkTrue((N - Gaussian(1.0, 1.0))$isapprox(Gaussian(24.000, 8.393),1e-3))
}
test_max_tuple = function(){
    checkEquals(c(0.1,0.05),max_tuple(c(0.,0.),c(0.1,0.05)))
    checkEquals(c(0.1,0.05),max_tuple(c(0.,0.05),c(0.1,0.0)))
}
test_list_diff = function(){
    old = list("0"=Gaussian(2.1,3.05))
    new = list("0"=Gaussian(2,3.0))
    checkEquals(list_diff(old,new),c(0.1,0.05))
}
test_1vs1 = function(){
    teams = list(c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300)) )
    result = c(1,0)
    g = Game(teams,result, 0.0)
    post = g$posteriors()
    checkTrue(post[[1]][[1]]$isapprox(Gaussian(20.7947792, 7.19448142)))
    checkTrue(post[[2]][[1]]$isapprox(Gaussian(29.2052207, 7.19448142)))
    
    teams = list("a"=c(Rating(29.,1.,25.0/6)),"b"= c(Rating(25.0,25.0/3,25.0/6)))
    g = Game(teams, c(1,0))
    teams <- g$posteriors()
    checkTrue(teams$a[[1]]$isapprox(Gaussian(28.89648,0.9966043)))
    checkTrue(teams$b[[1]]$isapprox(Gaussian(32.18921,6.062064)))
    
    teams = list("a"=c(Rating(1.139,0.531,1.0,0.2125)),"b"= c(Rating(15.568,0.51,1.0,0.2125)))
    g = Game(teams, c(1,0), 0.0)
    checkTrue(g$likelihoods[[1]][[1]]==Ninf)
    checkTrue(g$likelihoods[[2]][[1]]==Ninf)    
}
test_1vs1vs1 = function(){
    p = Game(list(c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300))), c(1,0,2))$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(25,6.238469796)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(31.31135822,6.6988186)))
    checkTrue(p[[3]][[1]]$isapprox(Gaussian(18.6886417,6.6988186)))
    
    p = Game(list(c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300))), c(1,0,2), 0.5)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(25,6.092561)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(33.37932,6.483576)))
    checkTrue(p[[3]][[1]]$isapprox(Gaussian(16.62068,6.483576)))
}
test_1vs1_draw = function(){
    p = Game(list(c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300))), c(0,0), 0.25)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(25,6.469481)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(25,6.469481)))
    
    p = Game(list(c(Rating(25.0,3,25.0/6,25.0/300)), c(Rating(29.0,2,25.0/6,25.0/300))), c(0,0), 0.25)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(25.736,2.709956)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(28.67289,1.916471)))
}
test_1vs1vs1_draw = function(){
    p = Game(list(c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300)), c(Rating(25.0,25.0/3,25.0/6,25.0/300))), c(0,0,0), 0.25)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(25,5.729069)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(25,5.707424)))
    
    p = Game(list(c(Rating(25.0,3,25.0/6,25.0/300)), c(Rating(25.0,3,25.0/6,25.0/300)), c(Rating(29.0,2,25.0/6,25.0/300))), c(0,0,0), 0.25)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(25.48851,2.638266)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(25.51067,2.628752)))
    checkTrue(p[[3]][[1]]$isapprox(Gaussian(28.55592,1.885689)))
}
test_NvsN_draw = function(){
    p = Game(list(c(Rating(15.0,1,25.0/6,25.0/300), Rating(15.0,1,25.0/6,25.0/300)), c(Rating(30,2,25.0/6,25.0/300))), c(0,0), 0.25)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(15,0.9916146)))
    checkTrue(p[[1]][[2]]$isapprox(Gaussian(15,0.9916146)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(30,1.932044)))
}
test_NvsNvsN_mixt = function(){
    teams = list(
    c(Rating(12.0,3,25.0/6,25.0/300), Rating(18.0,3,25.0/6,25.0/300)), c(Rating(30,3,25.0/6,25.0/300)),
    c(Rating(14.0,3,25.0/6,25.0/300), Rating(16.0,3,25.0/6,25.0/300)) )
    p = Game(teams, c(0,1,1), 0.25)$posteriors()
    checkTrue(p[[1]][[1]]$isapprox(Gaussian(13.05056,2.864398)))
    checkTrue(p[[1]][[2]]$isapprox(Gaussian(19.05056,2.864398)))
    checkTrue(p[[2]][[1]]$isapprox(Gaussian(29.29189,2.76353)))
    checkTrue(p[[3]][[1]]$isapprox(Gaussian(13.65756,2.813162)))
    checkTrue(p[[3]][[2]]$isapprox(Gaussian(15.65756,2.813162)))
}
test_game_evidence_1vs1 = function(){
    g = Game(list(c(Rating(25.0,1e-7,25.0/6,25.0/300)), c(Rating(25.0,1e-7,25.0/6,25.0/300))), c(0,0), 0.25)
    checkEquals(g$evidence,0.25)
    g = Game(list(c(Rating(25.0,1e-7,25.0/6,25.0/300)), c(Rating(25.0,1e-7,25.0/6,25.0/300))), c(0,1), 0.25)
    checkEquals(g$evidence,0.375)
}
test_game_evidence_1vs1vs1 = function(){
    teams = list(c(Rating(25.,1e-7,25.0/6,25.0/300)), c(Rating(25.,1e-7,25.0/6,25.0/300)), c(Rating(25.,1e-7,25.0/6,25.0/300)))
    
    g_abc = Game(teams, c(1,2,3), 0.)
    g_acb = Game(teams, c(1,3,2), 0.)
    g_bac = Game(teams, c(2,1,3), 0.)
    g_bca = Game(teams, c(3,1,2), 0.)
    g_cab = Game(teams, c(2,3,1), 0.)
    g_cba = Game(teams, c(3,2,1), 0.)
    proba = 0
    proba = proba + g_abc$evidence
    proba = proba + g_acb$evidence
    proba = proba + g_bac$evidence
    proba = proba + g_bca$evidence
    proba = proba + g_cab$evidence
    proba = proba + g_cba$evidence
    checkEquals(proba,1.49999991)
}
test_forget = function(){
    gamma = 0.15*25.0/3
    N = Gaussian(25,1e-7)
    checkEquals(N$forget(gamma,5)$sigma, sqrt(5*gamma^2))
    checkEquals(N$forget(gamma,1)$sigma, sqrt(1*gamma^2))
}
test_initialize_events = function(){
    composition = list(list(c("a"),c("b")),list(c("c"),c("d")),list(c("e"),c("f")))
    results = list(c(0,1),c(1,0),c(0,1))
    events = initialize_events(composition, results)
    checkTrue(events[[1]]$teams[[1]]$items[[1]]$name == "a")
}
test_initialize_skills = function(){
    composition = list(list(c("a"),c("b")),list(c("c"),c("d")),list(c("e"),c("f")))
    agents = list()
    for (a in c("a","b","c","d","e","f")){
        agents[[a]] = Agent(Rating(25,25/3,25/6,25/300), Ninf, -Inf)
    }
    
    skills = initialize_skills(composition, agents, 0)
    checkTrue(skills[[1]]$backward == Ninf)
}
test_one_event_each = function(){
    composition = list(list(c("a"),c("b")),list(c("c"),c("d")),list(c("e"),c("f")))
    results = list(c(0,1),c(1,0),c(0,1))
    agents = list()
    for (a in c("a","b","c","d","e","f")){
        agents[[a]] = Agent(Rating(25,25/3,25/6,25/300), Ninf, -Inf)
    }
    
    b = Batch(composition, results, 0, agents, Environment())
    p = b$posteriors()
    checkTrue(p$a == Gaussian(29.205,7.194))
    checkTrue(p$d == Gaussian(29.205,7.194))
    checkTrue(p$e == Gaussian(29.205,7.194))
    checkTrue(p$b == Gaussian(20.795,7.194))
    checkTrue(p$c == Gaussian(20.795,7.194))
}
test_batch_same_strength = function(){
    composition = list(list(c("a"),c("b")),list(c("a"),c("c")),list(c("b"),c("c")))
    results = list(c(0,1),c(1,0),c(0,1))
    agents = list()
    for (a in c("a","b","c")){
        agents[[a]] = Agent(Rating(25,25/3,25/6,25/300), Ninf, -Inf)
    }
    b = Batch(composition, results, 0, agents, Environment())
    p = b$posteriors()
    checkTrue(p$a == Gaussian(mu=24.961, sigma=6.299))
    checkTrue(p$b == Gaussian(mu=27.096, sigma=6.01))
    checkTrue(p$c == Gaussian(mu=24.89, sigma=5.866))
    b$convergence(1e-6,10)
    p = b$posteriors()
    checkTrue(p$a == Gaussian(25, 5.419))
    checkTrue(p$b == Gaussian(25, 5.419))
    checkTrue(p$c == Gaussian(25, 5.419))
    
    b$forward_prior_out("a")
    
}
test_history_init = function(){
    composition = list(list(c("aa"),c("b")),list(c("aa"),c("c")),list(c("b"),c("c")))
    results = list(c(0,1),c(1,0),c(0,1))
    priors = list()
    gamma = 0.15*25/3
    for (a in c("aa","b","c")){
        priors[[a]] = Rating(25,25/3,25/6,gamma)
    }
    h = History(composition, results, c(1,2,3), priors)
    p1 = h$batches[[1]]$posteriors()
    checkTrue(p1$aa$isapprox(Gaussian(29.205, 7.194),1e-3))
    checkTrue(p1$b$isapprox( Gaussian(20.795,7.194),1e-3))
    observed = h$batches[[2]]$skills[["aa"]]$forward$sigma
    checkEquals(observed,sqrt(gamma^2+ p1$aa$sigma^2))
    p2 = h$batches[[2]]$posteriors()
    checkTrue(p2$aa$isapprox(Gaussian(24.86, 6.374),1e-3))
    checkTrue(p2$c$isapprox(Gaussian(30.659, 6.922),1e-3))
}
test_one_batch_history = function(){
    composition = list(list(c("aj"),c("bj")),list(c("bj"),c("cj")),list(c("cj"),c("aj")))
    results = list(c(0,1),c(0,1),c(0,1))
    env = Environment(mu=25,sigma=25/3,beta=25/6,gamma=0.15*25/3)
    h = History(composition, results, times=c(0,0,0), env=env)
    p1 = h$batches[[1]]$posteriors()
    checkTrue(p1$aj$isapprox(Gaussian(22.904, 6.010),1e-3))
    checkTrue(p1$bj$isapprox(Gaussian(25.039, 6.299),1e-3))
    checkTrue(p1$cj$isapprox(Gaussian(25.11, 5.866),1e-3))
    step = h$convergence(verbose=T)
    p1 = h$batches[[1]]$posteriors()
    checkTrue(p1$aj$isapprox(Gaussian(25, 5.419),1e-3))
    checkTrue(p1$bj$isapprox(Gaussian(25, 5.419),1e-3))
    checkTrue(p1$cj$isapprox(Gaussian(25, 5.419),1e-3))

    env = Environment(mu=25,sigma=25/3,beta=25/6,gamma=25/300)
    h1 = History(composition=composition, results=results, times=c(1,2,3), env=env)
    p3 = h1$batches[[3]]$posteriors()
    checkTrue(p3$aj$isapprox(Gaussian(22.904, 6.011),1e-3))
    checkTrue(p3$cj$isapprox(Gaussian(25.11, 5.867),1e-3))
    step = h1$convergence(T)
    p3 = h1$batches[[3]]$posteriors()
    checkTrue(p3$aj$isapprox(Gaussian(24.999, 5.420),1e-3))
    checkTrue(p3$cj$isapprox(Gaussian(25.001, 5.420),1e-3))
}
test_trueSkill_Through_Time = function(){
    composition = list(list(c("a"),c("b")),list(c("a"),c("c")),list(c("b"),c("c")))
    results = list(c(0,1),c(1,0),c(0,1))
    env = Environment(mu=25,sigma=25/3,beta=25/6,gamma=25/300)
    h = History(composition=composition, results=results, times=c(), env=env)
    h$convergence()
    checkTrue(h$batches[[3]]$skills[["b"]]$elapsed==1)
    p1=h$batches[[1]]$posteriors(); p3=h$batches[[3]]$posteriors()
    checkTrue(p1$a$isapprox(Gaussian(25.000267, 5.4193816)))
    checkTrue(p1$b$isapprox(Gaussian(24.999465, 5.4194258)))
    checkTrue(p3$b$isapprox(Gaussian(25.00053219, 5.419696790)))
}
test_env_0_TTT = function(){
    composition = list(list(c("a"),c("b")),list(c("a"),c("c")),list(c("b"),c("c")))
    results = list(c(0,1),c(1,0),c(0,1))
    env = Environment(mu=0,sigma=6,beta=1,gamma=0.05)
    h = History(composition=composition, results=results, env=env)
    h$convergence()
    p1=h$batches[[1]]$posteriors(); p3=h$batches[[3]]$posteriors()
    checkTrue(p1$a$isapprox(Gaussian(0.001, 2.396),1e-3))
    checkTrue(p1$b$isapprox(Gaussian(-0.001,2.396),1e-3))
    checkTrue(p3$b$isapprox(Gaussian(0.001, 2.396),1e-3))
}
test_teams = function(){
    composition = list(list(c("a","b"),c("c","d")),list(c("e","f"),c("b","c")),list(c("a","d"),c("e","f")))
    results = list(c(0,1),c(1,0),c(0,1))
    env = Environment(mu=0,sigma=6,beta=1,gamma=0)
    h = History(composition=composition, results=results, env=env)
    h$convergence()
    p1=h$batches[[1]]$posteriors(); p2=h$batches[[2]]$posteriors()
    checkTrue(p1$a==p1$b)
    checkTrue(p1$c==p1$d)
    checkTrue(p2$f==p2$e)
    checkTrue(p1$a$isapprox(Gaussian(4.085,5.107),1e-3))
    checkTrue(p1$c$isapprox(Gaussian(-0.533,5.107),1e-3))
    p3=h$batches[[3]]$posteriors()
    checkTrue(p3$e$isapprox(Gaussian(-3.552,5.155),1e-3))
}
test_sigma_beta_0 = function(){
    composition = list(list(c("a","a_b","b"),c("c","c_d","d")),list(c("e","e_f","f"),c("b","b_c","c")),list(c("a","a_d","d"),c("e","e_f","f")))
    results = list(c(0,1),c(1,0),c(0,1))
    env = Environment(mu=0,sigma=6,beta=1,gamma=0)
    priors = list()
    for (a in c("a_b","c_d","e_f","b_c","a_d","e_f")){
        priors[[a]] = Rating(0,1e-7,0.0,0.2)
    }
    h = History(composition=composition, results=results, priors=priors, env=env)
    h$convergence()
    p1=h$batches[[1]]$posteriors()
    p2=h$batches[[2]]$posteriors()
    p3=h$batches[[3]]$posteriors()
    checkTrue(p1$a_b$isapprox(Gaussian(0,0),1e-4))
    checkTrue(p3$e_f$isapprox(Gaussian(-0.002,0.2),1e-4))
}
test_memory_size = function(){
    cat("Add a test for memory size please\n")
}
test_learning_curve = function(){
    cat("Add learning curve function\n")
}

source("TrueSkill.R")
if (!require("RUnit", quietly = TRUE)) {
  stop("Package Runit is not found.") 
}
