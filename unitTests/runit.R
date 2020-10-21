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

source("TrueSkill.R")
if (!require("RUnit", quietly = TRUE)) {
  stop("Package Runit is not found.") 
}
