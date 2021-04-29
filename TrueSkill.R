#library(microbenchmark)
library(hash)
library(compiler)

# Y really appreciate R, but: https://www.refsmmat.com/posts/2016-09-12-r-lists.html
# --
# look: https://bookdown.org/csgillespie/efficientR/programming.html
# --
# Environment: http://adv-r.had.co.nz/Environments.html
# library(pryr)

# --
# TODO: never to grow vectors. (vec <- rep(NA,10))
# TODO: hashable structure


BETA = 1.0
MU = 0.0
SIGMA = BETA * 6
GAMMA = BETA * 0.03
P_DRAW = 0.0
EPSILON = 1e-6
ITERATIONS = 30
PI = 1/(SIGMA^2)
TAU = MU*PI

sqrt2 = sqrt(2)

erfc <- function(x){
    z = abs(x)
    t = 1.0 / (1.0 + z / 2.0)
    a = -0.82215223 + t * 0.17087277 
    b =  1.48851587 + t * a
    c = -1.13520398 + t * b
    d =  0.27886807 + t * c
    e = -0.18628806 + t * d
    f =  0.09678418 + t * e
    g =  0.37409196 + t * f
    h =  1.00002368 + t * g
    r = t * exp(-z * z - 1.26551223 + t * h)
    if(x < 0){
        r = 2.0 - r
    }
    return(r)
}
erfc = cmpfun(erfc)

erfcinv <- function(y){#y=0.3
    if (y >= 2) return(-Inf)
    if (y < 0) stop("Argument must be nonnegative")
    if (y == 0) return(Inf)
    zero_point = (y < 1)
    if (! zero_point) y = 2 - y
    t = sqrt(-2 * log(y / 2.0))
    x = -0.70711 * ((2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t)
    for (h in c(0,1,2)){
        err = erfc(x) - y
        x = x + ( err / (1.12837916709551257 * exp(-(x^2)) - x * err))
    }
    if (zero_point){
        r = x
    }else{ 
        r = -x
    }
    return(r)
}
erfcinv = cmpfun(erfcinv)


mu_sigma <- function(tau_, pi_){
  if (pi_ > 0.0){
    sigma = sqrt(1/pi_)
    mu = tau_ / pi_
    return(c(mu,sigma))
  }
  if (pi_ + 1e-9 < 0.0){
    stop("Precision should be greater than 0")
  }
  return(c(0, Inf))
}
mu_sigma  = cmpfun(mu_sigma )


cdf = function(x, mu, sigma){
  z = -(x - mu) / (sigma * sqrt2)
  return(0.5*erfc(z) )
}
cdf = cmpfun(cdf)
pdf = function(x, mu, sigma){
  normalizer = (sqrt(2*pi) * sigma)^-1
  functional = exp( -((x - mu)^2) / (2*sigma ^2) ) 
  return(normalizer * functional)
}
pdf = cmpfun(pdf)
ppf = function(p, mu, sigma){
  return( mu - sigma * sqrt2 * erfcinv(2 * p))
}
ppf = cmpfun(ppf)
trunc = function(mu, sigma, margin, tie){
  if (!tie){
    alpha_ = (margin-mu)/sigma
    v = pdf(-alpha_,0,1) / cdf(-alpha_,0,1)
    w = v * (v + (-alpha_)) 
  }else{
    alpha_ = (-margin-mu)/sigma
    beta_ = ( margin-mu)/sigma
    v = (pdf(alpha_,0,1)-pdf(beta_,0,1))/(cdf(beta_,0,1)-cdf(alpha_,0,1))
    u = (alpha_*pdf(alpha_,0,1)-beta_*pdf(beta_,0,1))/(cdf(beta_,0,1)-cdf(alpha_,0,1))
    w =  - ( u - v**2 ) 
  }
  mu_trunc = mu + sigma * v
  sigma_trunc = sigma * sqrt(1-w)
  return(c(mu_trunc, sigma_trunc))
}
trunc = cmpfun(trunc)
approx = function(N, margin, tie){
  m_s = trunc(N@mu, N@sigma, margin, tie)
  N@mu = m_s[1]
  N@sigma = m_s[2]
  return(N)
}
approx = cmpfun(approx)
compute_margin = function(p_draw, sd){
  return(abs(ppf(0.5-p_draw/2, 0.0, sd )))
}
compute_margin  = cmpfun(compute_margin  )
max_tuple = function(t1, t2){
  return(c(max(t1[1],t2[1]), max(t1[2],t2[2])))
}
gr_tuple = function(tup, threshold){
    return( (tup[1] > threshold) | (tup[2] > threshold) )
}
sortperm = function(xs, decreasing = F){
  return(order(xs, decreasing = decreasing))
}

Gaussian <- function(mu=0, sigma=1){
    if(sigma>=0.0){
      return(new("Gaussian",mu=mu,sigma=sigma))
    }else{
      stop("Require: (sigma >= 0.0)")
    }
}

f_tau <- function(a){
    if (a@sigma > 0.0){return(a@mu*a@sigma^-2)
    }else{return(Inf)}
}
f_tau = cmpfun(f_tau)
f_pi <- function(a){
    if (a@sigma > 0.0){return(a@sigma^-2)
    }else{return(Inf)}
}
f_pi = cmpfun(f_pi)



gaussian <- setClass("Gaussian"
        , representation(mu = "numeric",sigma = "numeric"))
setMethod("show","Gaussian", function(object) { cat(paste0("Gaussian(mu=", round(object@mu,3), ", sigma=", round(object@sigma,3),")\n"))})
Pi <- function(a) 0
setGeneric("Pi")
setMethod("Pi", "Gaussian", function(a) if (a@sigma > 0.0){return(a@sigma^-2)
    }else{return(Inf)})
Tau <- function(a) 0
setGeneric("Tau")
setMethod("Tau", "Gaussian", function(a){
    if (a@sigma > 0.0){return(a@mu*a@sigma^-2)
    }else{return(Inf)}
  })
forget <- function(a,gamma,t) 0
setGeneric("forget")
setMethod("forget", c("Gaussian","numeric","numeric"), function(a,gamma,t){
    a@sigma = sqrt(a@sigma^2 + t*gamma^2)
    return(a)
  })
exclude <- function(N,M) 0
setGeneric("exclude")
setMethod("exclude", c("Gaussian","Gaussian"), 
  function(N,M){
    N@mu = N@mu-M@mu
    N@sigma = sqrt(N@sigma^2 - M@sigma^2)
    return(N)
  })
isapprox <- function(N, M, tol=1e-4) 0
setGeneric("isapprox")
setMethod("isapprox", c("Gaussian", "Gaussian", "numeric") , 
  function(N,M,tol=1e-4){
    return(abs(N@mu-M@mu)<tol & abs(N@sigma-M@sigma)<tol)
  })
delta <- function(N,M) 0
setGeneric("delta")
setMethod("delta", c("Gaussian", "Gaussian") , 
  function(N,M){
    return( c(abs(N@mu - M@mu) , abs(N@sigma - M@sigma)))
  })
setMethod("+", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    e1@mu = e1@mu + e2@mu
    e1@sigma = sqrt((e1@sigma^2) + (e2@sigma^2) )
    return(e1)
})
setMethod("-", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    e1@mu = e1@mu - e2@mu
    e1@sigma = sqrt((e1@sigma^2) + (e2@sigma^2) )
    return(e1)
})
setMethod("*", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    tau_ = f_tau(e1) + f_tau(e2); pi_ =  f_pi(e1) + f_pi(e2)
    if (pi_ > 0.0){
        e1@sigma = sqrt(1/pi_)
        e1@mu = tau_ / pi_
    }else{
        e1@mu  = 0
        e1@sigma= Inf
    }
    return(e1)
})
setMethod("/", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    tau_ = f_tau(e1) - f_tau(e2); pi_ =  f_pi(e1) - f_pi(e2)
    if (pi_ > 0.0){
        e1@sigma = sqrt(1/pi_)
        e1@mu = tau_ / pi_
    }else{
        e1@mu  = 0
        e1@sigma= Inf
    }
    return(e1)
})
setMethod("==", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    mu = e1@mu - e2@mu < 1e-3
    sigma = if (e2@sigma == Inf | e1@sigma == Inf) (e1@sigma==e2@sigma) else  (e1@sigma - e2@sigma < 1e-3)
    return(mu & sigma)
})
  

N01 = Gaussian(0,1)
N00 = Gaussian(0,0)
Ninf = Gaussian(0,Inf)
Nms = Gaussian(MU, SIGMA)

list_diff = function(old, new){
  step = c(0,0)
  for (a in names(old)){#a="0"
    step = max_tuple(step, delta(old[[a]], new[[a]]))
  }
  return(step)
}


Player <- function(prior=Nms, beta=BETA, gamma=GAMMA){
    return(new("Player", prior = prior, beta = beta, gamma = gamma))
}
player <- setClass("Player"
        , representation(prior = "Gaussian",beta = "numeric",gamma = "numeric"))
setMethod("show", "Player", 
  function(object){
    cat(paste0("Player(Gaussian(mu=", round(object@prior@mu,3), ", sigma=", round(object@prior@sigma,3),"), beta=",round(object@beta,3), ", gamma=",round(object@gamma,3)),")\n")
  })
performance <- function(a) 0
setGeneric("performance")
setMethod("performance", "Player", 
  function(a){
    N = a@prior
    N@sigma = sqrt(a@prior@sigma^2 + a@beta^2)
    return(N)
  })

# 
#   
# Player <- setRefClass("Player", 
#   fields = list(prior = "Gaussian"
#                ,beta = "numeric"
#                ,gamma = "numeric"
#                ,draw = "Gaussian")
# )
# Player$methods( 
#   initialize = function(prior=Gaussian(MU, SIGMA), beta=BETA, gamma=GAMMA){
#     prior <<- prior; beta <<- beta; gamma <<- gamma; draw <<- Ninf
#   },
#   show = function() { cat(paste0("Player(Gaussian(mu=", round(prior@mu,3), ", sigma=", round(prior@sigma,3),"), beta=",round(beta,3), ", gamma=",round(gamma,3)),")\n")},
#   performance = function(){
#     mu = prior@mu
#     sigma = sqrt(prior@sigma^2)
#     sigma = sigma  + beta^2
#     return(Gaussian(mu, sigma))
#   }
# )

team_messages <- function(prior=Ninf, likelihood_lose=Ninf, likelihood_win=Ninf){
    return(new("team_messages", prior=prior, likelihood_lose = likelihood_lose,  likelihood_win = likelihood_win))
   }
Team_messages <- setClass("team_messages"
  , representation(prior = "Gaussian",
    likelihood_lose = "Gaussian",
    likelihood_win = "Gaussian")
)
posterior_win  <- function(object) 0
setGeneric("posterior_win")
setMethod("posterior_win", "team_messages", function(object){
  return(object@prior*object@likelihood_lose)
})
posterior_lose  <- function(object) 0
setGeneric("posterior_lose")
setMethod("posterior_lose", "team_messages", function(object){
  return(object@prior*object@likelihood_win)
})
likelihood <- function(object) 0
setGeneric("likelihood")
setMethod("likelihood", "team_messages", function(object){
  return(object@likelihood_win*object@likelihood_lose)
  })


diff_messages <- function(prior =Ninf, likelihood=Ninf){
    return(new("diff_messages", prior=prior, likelihood=likelihood))
}
Diff_messages <- setClass("diff_messages",
  representation( prior = "Gaussian",likelihood = "Gaussian")
  )

Game <- function(teams, result = vector(), p_draw=P_DRAW){
if ((length(result)>0) & (length(teams) != length(result))) stop("(length(results)>0) & (length(teams) != length(result))")
    if ((0.0 > p_draw) | (1.0 <= p_draw)) stop("0.0 <= p_draw < 1.0")
    if ((p_draw==0.0) & (length(result)>0) & (length(unique(result))!=length(result))) stop("(p_draw=0.0) & (length(result)>0) & (length(unique(result))!=length(result))")
    
    l_e = compute_likelihoods(teams, result, p_draw)
    g = new("Game", teams= teams, result = result, p_draw = p_draw, likelihoods = l_e$likelihoods, evidence = l_e$evidence)
    return(g)
}
game <- setClass("Game",
  representation(
    teams = "list",
    result = "vector",
    p_draw = "numeric",
    likelihoods = "list",
    evidence = "numeric")
)
# setMethod("size", "Game", function(object){
#     res = rep(NA,length(object@teams))
#     for (t in seq(length(object@teams))){ res[t] =length(team[t]) }
#     return(res)
#   })
graphical_model <- function(teams, result, p_draw){
    r = if (length(result)>0) result else seq(length(teams)-1,0)
    o = sortperm(r, decreasing = T)
    t = vector('list', length(teams))
    d = vector('list', length(teams)-1)
    tie = rep(NA, length(d)); margin =  rep(NA, length(d))
    for (e in seq(length(teams))){#e=1
      team_perf = N00
      for (a in teams[[o[e]]]){ team_perf  = team_perf + performance(a)}
      t[[e]] = team_messages(team_perf) }
    for (e in seq(length(teams)-1)){ 
      d[[e]] = diff_messages(t[[e]]@prior - t[[e+1]]@prior) }
    for (e in seq(length(d))){
     tie[e] = r[o[e]]==r[o[e+1]]}
    for (e in seq(length(d))){
      if (p_draw == 0.0){ 
        margin[e] = 0.0}
      else{ 
        betas = 0
        for (a in teams[[o[e]]]){betas = betas + a@beta^2}
        for (a in teams[[o[e+1]]]){betas = betas + a@beta^2}
        margin[e] = compute_margin(p_draw, sqrt(betas) ) 
      }
    }
    evidence = 1
    for (e in seq(length(d))){
        mu = d[[e]]@prior@mu; sigma = d[[e]]@prior@sigma
        evidence = evidence * (if (tie[e]) (cdf(margin[e],mu,sigma)-cdf(-margin[e],mu,sigma)) else 1-cdf(margin[e],mu,sigma))
    }
    return(list("o"=o, "t"=t, "d"=d, "tie"=tie, "margin"=margin, "evidence" = evidence))
}
graphical_model = cmpfun(graphical_model )
likelihood_analitico <- function(teams,result,p_draw,gr){
    d = gr$d[[1]]@prior
    mu_sigma_trunc =  trunc(d@mu, d@sigma, gr$margin[1], gr$tie[1])
    mu_trunc = mu_sigma_trunc[1]; sigma_trunc = mu_sigma_trunc[2]
    if (d@sigma==sigma_trunc){
      delta_div = 0.0
      theta_div_pow2 = Inf
    }else{
      delta_div = (d@sigma^2*mu_trunc - sigma_trunc^2*d@mu)/(d@sigma^2-sigma_trunc^2)
      theta_div_pow2 = (sigma_trunc^2*d@sigma^2)/(d@sigma^2 - sigma_trunc^2)
    }
    res = vector('list', length(teams))
    for (i in seq(length(teams))){#i=1
      team = vector('list', length(teams[[i]]))
      for (j in seq(length(teams[[i]]))){#j=1
        N = d
        N@mu = if (d@sigma == sigma_trunc) 0.0 else teams[[i]][[j]]@prior@mu + ( delta_div - d@mu)*(-1)^(gr$o[i]==2)
        N@sigma = sqrt(theta_div_pow2 + d@sigma^2 - teams[[i]][[j]]@prior@sigma^2)
        team[[j]] = N
      }
      res[[i]] = team
    }
    return(res)
}
likelihood_analitico = cmpfun(likelihood_analitico )
likelihood_teams <- function(teams,result,p_draw,gr){
    o= gr$o; d = gr$d; t = gr$t; margin = gr$margin; tie = gr$tie 
    step = c(Inf,Inf); i = 0
    while (gr_tuple(step,1e-3) & i < 10){
      for (e in seq(length(d)-1)){
        d[[e]]@prior = posterior_win(t[[e]]) - posterior_lose(t[[e+1]])
        d[[e]]@likelihood = approx(d[[e]]@prior,margin[[e]],tie[[e]])/d[[e]]@prior
        likelihood_lose = posterior_win(t[[e]]) - d[[e]]@likelihood
        step = max_tuple(step,delta(t[[e+1]]@likelihood_lose,likelihood_lose))
        t[[e+1]]@likelihood_lose = likelihood_lose
      }
      for (e in seq(length(d),2,-1)){
        d[[e]]@prior = posterior_win(t[[e]]) - posterior_lose(t[[e+1]])
        d[[e]]@likelihood = approx(d[[e]]@prior,margin[[e]],tie[[e]])/d[[e]]@prior
        likelihood_win = posterior_lose(t[[e+1]]) + d[[e]]@likelihood
        step = max_tuple(step,delta(t[[e]]@likelihood_win,likelihood_win))
        t[[e]]@likelihood_win = likelihood_win
      }
      i = i + 1
    }
    if (length(d)==1){
      d[[1]]@prior = posterior_win(t[[1]]) - posterior_lose(t[[2]])
      d[[1]]@likelihood = approx(d[[1]]@prior,margin[[1]],tie[[1]])/d[[1]]@prior
    }
    t[[1]]@likelihood_win = posterior_lose(t[[2]]) + d[[1]]@likelihood
    t[[length(t)]]@likelihood_lose = posterior_win(t[[length(t)-1]]) - d[[length(d)]]@likelihood
    res = vector('list', length(t))
    for (e in seq(length(t))){ res[[e]] = likelihood(t[[o[e]]])}
    return(res)
}
likelihood_teams = cmpfun(likelihood_teams )
compute_likelihoods <-  function(teams,result,p_draw){
    gr = graphical_model(teams,result,p_draw)
    if (length(teams)>2){
      lhoods = vector('list', length(teams))
      m_t_ft = likelihood_teams(teams,result,p_draw,gr)
      for (e in seq(length(teams))){#e=2
        lhoods[[e]]  = vector('list', length(teams[[e]]))
        for (i in seq(length(teams[[e]]))){#i=1
          team_perf = N00
          for (a in teams[[e]]){ team_perf  = team_perf + performance(a)}
          ex = exclude(team_perf,teams[[e]][[i]]@prior)
          lhoods[[e]][[i]] = m_t_ft[[e]] - ex
        }
      }
      return(list("likelihoods"=lhoods, "evidence"=gr$evidence))
    }else{
      return(list("likelihoods"=likelihood_analitico(teams,result,p_draw,gr), "evidence"=gr$evidence))
    }
}
compute_likelihoods = cmpfun(compute_likelihoods)
posteriors <- function(g) 0
setGeneric("posteriors")
setMethod("posteriors", "Game", function(g){
    res = vector('list', length(g@teams))
    for (e in seq(length(g@teams))){
      post = vector('list', length(g@teams[[e]]))
      for (i in seq(length(g@teams[[e]]))){
        post[[i]] = g@teams[[e]][[i]]@prior * g@likelihoods[[e]][[i]]
      }
      if (is.null(names(g@teams))){
        res[[e]] = post
      }else{
        res[[names(g@teams)[e]]] = post
      }
    }
    return(res)
  }
)

# 
# 
#   
#   
#   
#   
# Game$methods(
#   size = function(){
#     res = rep(NA,length(teams))
#     for (t in seq(length(teams))){ res[t] =length(team[t]) }
#     return(res)
#   },
#   graphical_model = function(){
#     r = if (length(result)>0) result else seq(length(teams)-1,0)
#     o = sortperm(r, decreasing = T)
#     t = vector('list', length(teams))
#     d = vector('list', length(teams)-1)
#     tie = rep(NA, length(d)); margin =  rep(NA, length(d))
#     for (e in seq(length(teams))){#e=1
#       team_perf = Gaussian(0,0)
#       for (a in teams[[o[e]]]){ team_perf  = team_perf + performance(a)}
#       t[[e]] = team_messages(team_perf) }
#     for (e in seq(length(teams)-1)){ 
#       d[[e]] = diff_messages(t[[e]]@prior - t[[e+1]]@prior) }
#     for (e in seq(length(d))){
#      tie[e] = r[o[e]]==r[o[e+1]]}
#     for (e in seq(length(d))){
#       if (p_draw == 0.0){ 
#         margin[e] = 0.0}
#       else{ 
#         betas = 0
#         for (a in teams[[o[e]]]){betas = betas + a@beta^2}
#         for (a in teams[[o[e+1]]]){betas = betas + a@beta^2}
#         margin[e] = compute_margin(p_draw, sqrt(betas) ) 
#       }
#     }
#     evidence <<- 1
#     for (e in seq(length(d))){
#         mu = d[[e]]@prior@mu; sigma = d[[e]]@prior@sigma
#         evidence <<- evidence * (if (tie[e]) (cdf(margin[e],mu,sigma)-cdf(-margin[e],mu,sigma)) else 1-cdf(margin[e],mu,sigma))
#     }
#     return(list("o"=o, "t"=t, "d"=d, "tie"=tie, "margin"=margin))
#   },
#   likelihood_analitico = function(){
#     gr = graphical_model()
#     d = gr$d[[1]]@prior
#     mu_sigma_trunc =  trunc(d@mu, d@sigma, gr$margin[1], gr$tie[1])
#     mu_trunc = mu_sigma_trunc[1]; sigma_trunc = mu_sigma_trunc[2]
#     if (d@sigma==sigma_trunc){
#       delta_div = 0.0
#       theta_div_pow2 = Inf
#     }else{
#       delta_div = (d@sigma^2*mu_trunc - sigma_trunc^2*d@mu)/(d@sigma^2-sigma_trunc^2)
#       theta_div_pow2 = (sigma_trunc^2*d@sigma^2)/(d@sigma^2 - sigma_trunc^2)
#     }
#     res = vector('list', length(teams))
#     for (i in seq(length(teams))){#i=1
#       team = vector('list', length(teams[[i]]))
#       for (j in seq(length(teams[[i]]))){#j=1
#         mu = if (d@sigma == sigma_trunc) 0.0 else teams[[i]][[j]]@prior@mu + ( delta_div - d@mu)*(-1)^(gr$o[i]==2)
#         sigma_analitico = sqrt(theta_div_pow2 + d@sigma^2 - teams[[i]][[j]]@prior@sigma^2)
#         team[[j]] = Gaussian(mu,sigma_analitico)
#       }
#       res[[i]] = team
#     }
#     return(res)
#   },
#   likelihood_teams = function(){
#     gr = graphical_model()
#     o= gr$o; d = gr$d; t = gr$t; margin = gr$margin; tie = gr$tie 
#     step = c(Inf,Inf); i = 0
#     while (gr_tuple(step,1e-6) & i < 10){
#       for (e in seq(length(d)-1)){
#         d[[e]]@prior = posterior_win(t[[e]]) - posterior_lose(t[[e+1]])
#         d[[e]]@likelihood = approx(d[[e]]@prior,margin[[e]],tie[[e]])/d[[e]]@prior
#         likelihood_lose = posterior_win(t[[e]]) - d[[e]]@likelihood
#         step = max_tuple(step,delta(t[[e+1]]@likelihood_lose,likelihood_lose))
#         t[[e+1]]@likelihood_lose = likelihood_lose
#       }
#       for (e in seq(length(d),2,-1)){
#         d[[e]]@prior = posterior_win(t[[e]]) - posterior_lose(t[[e+1]])
#         d[[e]]@likelihood = approx(d[[e]]@prior,margin[[e]],tie[[e]])/d[[e]]@prior
#         likelihood_win = posterior_lose(t[[e+1]]) + d[[e]]@likelihood
#         step = max_tuple(step,delta(t[[e]]@likelihood_win,likelihood_win))
#         t[[e]]@likelihood_win = likelihood_win
#       }
#       i = i + 1
#     }
#     if (length(d)==1){
#       d[[1]]@prior = posterior_win(t[[1]]) - posterior_lose(t[[2]])
#       d[[1]]@likelihood = approx(d[[1]]@prior,margin[[1]],tie[[1]])/d[[1]]@prior
#     }
#     t[[1]]@likelihood_win = posterior_lose(t[[2]]) + d[[1]]@likelihood
#     t[[length(t)]]@likelihood_lose = posterior_win(t[[length(t)-1]]) - d[[length(d)]]@likelihood
#     res = vector('list', length(t))
#     for (e in seq(length(t))){ res[[e]] = likelihood(t[[o[e]]])}
#     return(res)
#   },
#   compute_likelihoods = function(){
#     if (length(teams)>2){
#       m_t_ft = likelihood_teams()
#       for (e in seq(length(teams))){#e=1
#         res = vector('list', length(teams[[e]]))
#         for (i in seq(length(teams[[e]]))){#i=1
#           team_perf = Gaussian(0,0)
#           for (a in teams[[e]]){ team_perf  = team_perf + performance(a)}
#           ex = exclude(team_perf,teams[[e]][[i]]@prior)
#           res[[i]] = m_t_ft[[e]] - ex
#         }
#         likelihoods[[e]] <<- res
#       }
#     }else{
#       likelihoods <<- likelihood_analitico()
#     }
#   },
#   posteriors = function(){
#     res = vector('list', length(teams))
#     for (e in seq(length(teams))){
#       post = vector('list', length(teams[[e]]))
#       for (i in seq(length(teams[[e]]))){
#         post[[i]] = likelihoods[[e]][[i]] * teams[[e]][[i]]@prior
#       }
#       if (is.null(names(teams))){
#         res[[e]] = post
#       }else{
#         res[[names(teams)[e]]] = post
#       }
#     }
#     return(res)
#   }
# )

Skill <- setRefClass("Skill",
  fields = list(
    forward = "Gaussian",
    backward = "Gaussian",
    likelihood = "Gaussian",
    elapsed = "numeric")
)
Skill$methods(
  initialize = function(forward=Ninf, backward=Ninf, likelihood=Ninf, elapsed=0){
    forward <<- forward; backward <<- backward
    likelihood <<- likelihood; elapsed <<- elapsed
  },
  posterior = function(){
    return(forward*likelihood*backward)
  },
  posterior_back = function(){
    return(forward*likelihood)
  },
  posterior_for = function(){
    return(likelihood*backward)
  }
)

Agent <- setRefClass("Agent",
  fields = list(
    player = "Player",
    message = "Gaussian",
    last_time = "numeric")
)
Agent$methods(
  initialize = function(player, message, last_time){
    player <<- player; message <<- message; last_time <<- last_time
  },
  receive = function(elapsed){
    if (!(message==Ninf)){
      res = forget(message,player@gamma, elapsed)
    }else{
      res = player@prior
    }
    return(res)
  }
)

Item <- setRefClass("Item",
  fields = list(name = "character", likelihood = "Gaussian")
)
Item$methods(  
  initialize = function(name, likelihood){
    name <<- name; likelihood <<- likelihood }
)

Team <- setRefClass("Team",
  fields = list(items = "vector", output = "numeric")
)
Team$methods(
  initialize = function(items, output){
    items <<- items; output <<- output }
)

Event <- setRefClass("Event",
  fields = list(teams = "vector", evidence = "numeric")
)
 Event$methods(
  initialize = function(teams, evidence){
    teams <<- teams
    evidence <<- evidence
  },
  #show = function(){
  #  print(paste("Event(",names(),",",result(),")"))
  #},
  names = function(){
    res = vector('list', length(teams))
    for (t in seq(length(teams))){
      vec = rep(NA, length(teams[[t]]))
      for (i in seq(length(teams[[t]]))){
        vec[i] = teams[[t]]$items[[i]]$name
      }
      res[[t]] = vec
    }
    return(res)
  },
  result = function(){
    res = c()
    for (team in teams){
      res = c(res, team$output)
    }
    return(res)
  }
)

compute_elapsed = function(last_time, actual_time){
  return(if (last_time == -Inf) 0 else (if (last_time == Inf) 1 else (actual_time - last_time)))
}

list_unique = function(xss){
  return(unique(unlist(xss)))
}

initialize_events = function(composition, results){
    events_ = vector('list',length(composition))
    for (e in seq(length(composition))){
      teams_ =  vector('list',length(composition[[e]]))
      for (t in seq(length(composition[[e]]))){
        items_ = vector('list',length(composition[[e]][[t]]))
        for (a in seq(length(composition[[e]][[t]]))){
          items_[[a]] = Item(composition[[e]][[t]][[a]], Ninf)
        }
        teams_[[t]] = Team(items_, if (length(results)>0) results[[e]][[t]] else length(composition[[e]])-t)
      }
      events_[[e]] = Event(teams_, 0)
    }
    return(events_)
}
initialize_skills = function(composition,agents,time){
    this_agents = list_unique(composition)
    skills_ = hash()
    for (a in this_agents){
      elapsed = compute_elapsed(agents[[a]]$last_time, time) 
      skills_[[a]] = Skill(agents[[a]]$receive(elapsed),Ninf,Ninf,elapsed)
    }
    return(skills_)
}

Batch <- setRefClass("Batch",
  fields = list(
    time = "numeric",
    events = "vector",
    skills = "hash",
    agents = "hash",
    p_draw = "numeric"
    )
)
Batch$methods(
  initialize = function(composition, results = list() ,time = 0, agents = hash(), p_draw=P_DRAW){
    if ((length(results)>0) & (length(composition) != length(results))) stop("(length(results)>0) & (length(composition) != length(results))")
      
    skills <<- initialize_skills(composition, agents, time)
    events <<- initialize_events(composition, results)
    time <<- time
    agents <<- agents
    p_draw <<- p_draw
    iteration()
  },
  show = function(){
    cat(paste0("Batch(time=",time,", events=",length(events),", skills=", length(skills),")\n"))
  },
  posterior = function(a){
    return(skills[[a]]$posterior() )
  },
  posteriors = function(){
    res = hash()
    for (a in names(skills)){
      res[[a]] = skills[[a]]$posterior() 
    }
    return(res)
  },
  within_prior = function(item){
    r = agents[[item$name]]$player
    r@prior = skills[[item$name]]$posterior()/item$likelihood
    return(r)
  },
  within_priors = function(event){
    res = list()
    for (t in seq(length(events[[event]]$teams))){
      vec = c()
      for (item in events[[event]]$teams[[t]]$items){
        vec = c(vec,within_prior(item))
      }
      res[[t]] = vec
    }
    return(res)
  },
  iteration = function(from_ =1){
    for (e in seq(from_,length(events))){
      teams = within_priors(e)
      result = events[[e]]$result()
      g = Game(teams, result, p_draw)
      t = 1
      for (team in events[[e]]$teams){
        i = 1
        for (item in team$items){
          skills[[item$name]]$likelihood <- skills[[item$name]]$likelihood/item$likelihood * g@likelihoods[[t]][[i]]
          item$likelihood = g@likelihoods[[t]][[i]]
          i = i + 1
        }
        t = t + 1
      }
      events[[e]]$evidence <- g@evidence
    }
  },
  convergence = function(epsilon,iterations){
    step = c(Inf,Inf); i = 0
    while (gr_tuple(step,epsilon) & i < iterations){
      old = posteriors()
      iteration()
      step = list_diff(old,posteriors())
      i = i + 1
    }
    return(step)
  },
  forward_prior_out = function(name){
    return(skills[[name]]$posterior_back())
  },
  backward_prior_out = function(name){
    return(forget(skills[[name]]$posterior_for(),agents[[name]]$player@gamma,skills[[name]]$elapsed))
  },
  new_backward_info = function(){
    for (a in names(skills)){
      skills[[a]]$backward <- agents[[a]]$message
    }
    return(iteration())
  },
  new_forward_info = function(){
    for (a in names(skills)){
      skills[[a]]$forward <- agents[[a]]$receive(skills[[a]]$elapsed)
    }
    return(iteration())
  }
)


History = setRefClass("History",
  fields = list(
    size = "numeric",
    batches = "vector",
    agents = "hash",
    time = "logical",
    mu = "numeric",
    sigma = "numeric",
    beta = "numeric",
    gamma = "numeric",
    p_draw = "numeric",
    h_epsilon = "numeric",
    h_iterations = "numeric"
    )
)
History$methods(
  initialize = function(composition, results=list(), times=c(), priors=hash(), mu=MU, sigma=SIGMA, beta=BETA, gamma=GAMMA, p_draw=P_DRAW, epsilon=EPSILON, iterations=ITERATIONS){
    if ((length(results)>0) & (length(composition) != length(results))){ stop("(length(results)>0) & (length(composition) != length(results))")}
    if (length(times) > 0 & (length(composition) != length(times))){ stop("length(times) error")}
    
    N_default = Gaussian(mu, sigma)
    this_agents = list_unique(composition)
    agents_ = hash()
    for (a in this_agents ){
        agents_[[a]] = Agent(if (a %in% names(priors)) priors[[a]] else Player(N_default, beta, gamma), Ninf, -Inf)
    }
    
     size <<- length(composition)
     agents <<- agents_
     time <<- length(times) > 0
     mu <<- mu; sigma <<- sigma; beta <<- beta; gamma <<- gamma; p_draw <<- p_draw; h_epsilon <<- epsilon; h_iterations <<- iterations
     trueskill(composition, results, times)
  },
  show = function(){
    cat(paste0("History(Events=",size, ", Batches=", length(batches), ", Agents=", length(agents),")\n"))
  },
  trueskill = function(composition, results, times){
    o = if (time) sortperm(times) else seq(size)
    i = 1
    while (i <= size){
      j = i; t = if (!time) i else times[o[i]]
      while (((time) & (j < size)) && (times[o[j+1]] == t)){j=j+1}
      if (length(results)>0){
        b = Batch(composition[o[seq(i,j)]], results[o[seq(i,j)]], t, agents, p_draw)
      }else{
        b = Batch(composition[o[seq(i,j)]], list(), t, agents, p_draw)  
      }
      batches <<- if (is.null(batches)) c(b) else c(batches, b)
      for (a in names(b$skills)){
        agents[[a]]$last_time <- if (!time) Inf else t
        agents[[a]]$message <- b$forward_prior_out(a)
      }
      i = j + 1
    }    
  },
  iteration = function(){
    step = c(0,0)
    
    if (length(batches)==1){
      old = batches[[1]]$posteriors()
      batches[[1]]$iteration()
      step = max_tuple(step,list_diff(old, batches[[1]]$posteriors()))
    
    }else{
    
    #clean(agents)
    for (a in names(agents)){agents[[a]]$message <<- Ninf}
    for (j in seq(length(batches)-1,1,-1)){
      for (a in names(batches[[j+1]]$skills)){
        agents[[a]]$message <- batches[[j+1]]$backward_prior_out(a)
      }
      old = batches[[j]]$posteriors()
      batches[[j]]$new_backward_info()
      step = max_tuple(step,list_diff(old, batches[[j]]$posteriors()))
    }
    #clean(agents)
    for (a in names(agents)){agents[[a]]$message <- Ninf}
    for (j in seq(2,length(batches))){
      for (a in names(batches[[j-1]]$skills)){
        agents[[a]]$message <- batches[[j-1]]$forward_prior_out(a)
      }
      old = batches[[j]]$posteriors()
      batches[[j]]$new_forward_info()
      step = max_tuple(step,list_diff(old, batches[[j]]$posteriors()))
    }
    
    } #end ELSE
    return(step)
  },
  convergence = function(epsilon = NA, iterations = NA, verbose=TRUE){
    if(is.na(epsilon)){epsilon = h_epsilon}
    if(is.na(iterations)){iterations = h_iterations}
    step = c(Inf, Inf); i = 1
    while (gr_tuple(step,epsilon) & i <= iterations){
      if (verbose){cat(paste0("Iteration = ", i))}
      step = iteration()
      i = i + 1
      if (verbose){cat(paste0(" step = (", step[1], ", ", step[2], ")\n"))}
    }
    if (verbose){cat("End\n")}
    return(step)
  },
  learning_curves = function(){
    res = hash()
    for (b in batches){
      for (a in names(b$skills)){
        t_p = list(t=b$time, N=b$posterior(a))
        if (has.key(a, res)){
          i = length(res[[a]])
          res[[a]][[i+1]] = t_p
        }else{
          res[[a]][[1]] = t_p
        }
      }
    }
    return(res)
  }
)

lc_print <- function(lc.a){
  res = "["
  for (i in seq(length(lc.a))){
    if (i == 1){
      res = paste0(res,"(",lc.a[[i]][[1]],", Gaussian(mu=",round(lc.a[[i]][[2]]@mu,3),", sigma=",round(lc.a[[i]][[2]]@sigma,3),"))")
    }else{
      res = paste0(res,", (",lc.a[[i]][[1]],", Gaussian(mu=",round(lc.a[[i]][[2]]@mu,3),", sigma=",round(lc.a[[i]][[2]]@sigma,3),"))")
    }
  }
  res = paste0(res,"]\n")
  cat(res)
}


# 
# teams = list(ta = c(Rating()), tb = c(Rating()))
# result = c(0,1)
# g = Game(teams, result)
# g$size()
# microbenchmark(Rating(0.001,1.0), times=5, unit="s")
# microbenchmark(N1/N2, times=5, unit="s")
# 
# 
