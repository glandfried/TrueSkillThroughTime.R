library(hash)
library(methods)
library(stats)

BETA = 1.0
MU = 0.0
SIGMA = BETA * 6
GAMMA = BETA * 0.03
P_DRAW = 0.0
EPSILON = 1e-6
ITERATIONS = 30

sqrt2 = sqrt(2)

cdf = function(x, mu, sigma){
  z = (x - mu) / (sigma)
  return(pnorm(z))
}
pdf = function(x, mu, sigma){
  return(dnorm(x,mu,sigma))
}
ppf = function(p, mu, sigma){
  return(qnorm(p, mu, sigma))
}

compute_margin = function(p_draw, sd){
  return(abs(ppf(0.5-p_draw/2, 0.0, sd )))
}

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

approx = function(N, margin, tie){
  m_s = trunc(N@mu, N@sigma, margin, tie)
  N@mu = m_s[1]
  N@sigma = m_s[2]
  return(N)
}
max_tuple = function(t1, t2){
  return(c(max(t1[1],t2[1]), max(t1[2],t2[2])))
}
gr_tuple = function(tup, threshold){
    return( (tup[1] > threshold) | (tup[2] > threshold) )
}
sortperm = function(xs, decreasing = F){
  return(order(xs, decreasing = decreasing))
}

#' @title Gaussian
#'
#' @description Gaussian class
#'
#' @param mu A number, the mean of the Gaussian distribution.
#' @param sigma A number, the standar deviation of the Gaussian distribution.
#' @param N A Gaussian object
#' @param M A Gaussian object
#' @param gamma The dynamic factor, the dynamic uncertainty
#' @param tol The tolerance threshold for comparitions
#' @param e1 A Gaussian object
#' @param e2 A Gaussian object
#' @param a A Gaussian object
#' @param t The elapsed time 
#'
#' @return Gaussian object
#'
#' @examples
#' N01 = Gaussian(0,1); N12 = Gaussian(mu = 1, sigma = 2)
#' N06 = Gaussian(); Ninf = Gaussian(0,Inf)
#' N01 * Ninf == N01
#' N01 * N12
#' N01 / N12
#' N01 + N12
#' N01 - N12
#' Pi(N12) == 1/(N12@sigma^2)
#' Tau(N12) == N12@mu/(N12@sigma^2)
#' Nnew = forget(N = N01, gamma = 0.01, t = 100)
#' isapprox(Nnew, Gaussian(N01@mu,sqrt(N01@sigma^2+100*(0.01^2))), tol=1e-6)
#'
#' @name Gaussian
#' @import methods
#' @import stats
#' @export
Gaussian <- function(mu=0, sigma=1){
    if(sigma>=0.0){
      return(new("Gaussian",mu=mu,sigma=sigma))
    }else{
      stop("Require: (sigma >= 0.0)")
    }
}
f_tau <- function(e1){
    if (e1@sigma > 0.0){return(e1@mu*e1@sigma^-2)
    }else{return(Inf)}
}
f_pi <- function(e1){
    if (e1@sigma > 0.0){return(e1@sigma^-2)
    }else{return(Inf)}
}
gaussian <- setClass("Gaussian"
        , representation(mu = "numeric",sigma = "numeric"))
setMethod("show","Gaussian", function(object) { cat(paste0("Gaussian(mu=", round(object@mu,3), ", sigma=", round(object@sigma,3),")\n"))})
#' @rdname Gaussian
#' @export
Pi <- function(N) 0
setGeneric("Pi")
#' @rdname Gaussian
#' @export
setMethod("Pi", "Gaussian", function(N) if (N@sigma > 0.0){return(N@sigma^-2)
    }else{return(Inf)})
#' @rdname Gaussian
#' @export
Tau <- function(N) 0
setGeneric("Tau")
#' @rdname Gaussian
#' @export
setMethod("Tau", "Gaussian", function(N){
    if (N@sigma > 0.0){return(N@mu*N@sigma^-2)
    }else{return(Inf)}
  })
#' @rdname Gaussian
#' @export
forget <- function(N,gamma,t) 0
setGeneric("forget")
#' @rdname Gaussian
#' @export
setMethod("forget", c("Gaussian","numeric","numeric"), function(N,gamma,t){
    N@sigma = sqrt(N@sigma^2 + t*gamma^2)
    return(N)
  })
exclude <- function(N,M) 0
setGeneric("exclude")
setMethod("exclude", c("Gaussian","Gaussian"), 
  function(N,M){
    N@mu = N@mu-M@mu
    N@sigma = sqrt(N@sigma^2 - M@sigma^2)
    return(N)
  })
#' @rdname Gaussian
#' @export
isapprox <- function(N, M, tol=1e-4) 0
setGeneric("isapprox")
#' @rdname Gaussian
#' @export
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
#' @rdname Gaussian
#' @export
setMethod("+", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    e1@mu = e1@mu + e2@mu
    e1@sigma = sqrt((e1@sigma^2) + (e2@sigma^2) )
    return(e1)
})
#' @rdname Gaussian
#' @export
setMethod("-", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    e1@mu = e1@mu - e2@mu
    e1@sigma = sqrt((e1@sigma^2) + (e2@sigma^2) )
    return(e1)
})
#' @rdname Gaussian
#' @export
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
#' @rdname Gaussian
#' @export
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
#' @rdname Gaussian
#' @export
setMethod("==", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    mu = abs(e1@mu - e2@mu) < 1e-3
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

#' @title Player
#'
#' @description Player class
#'
#' @param prior A Gaussian object, the prior belief distribution of the skills. The 
#' default value is: \code{Nms = Gaussian(mu = 0, sigma = 6)}.
#' @param beta A number, the standard deviation of the performance. The default 
#' value is: \code{BETA = 1}. The parameter \code{beta} acts as the scale of the 
#' estimates. A real difference of one \code{beta} between two skills is equivalent 
#' to 76\% probability of winning.
#' @param gamma A number, the amount of uncertainty (standar deviation) added to 
#' the estimates at each time step. The default value is: \code{GAMMA = 0.03}.
#' @param a A Player object
#'
#' @return Player object
#'
#' @examples
#' a1 = Player(prior = Gaussian(0,6), beta = 1, gamma = 0.03); 
#' a2 = Player()
#' a1@gamma == a2@gamma 
#' N = performance(a1) 
#' N@mu == a1@prior@mu
#' N@sigma == sqrt(a1@prior@sigma^2 + a1@beta^2)
#'
#' @name Player
#' @export
Player <- function(prior=Nms, beta=BETA, gamma=GAMMA){
    return(new("Player", prior = prior, beta = beta, gamma = gamma))
}
player <- setClass("Player"
        , representation(prior = "Gaussian",beta = "numeric",gamma = "numeric"))
setMethod("show", "Player", 
  function(object){
    cat(paste0("Player(Gaussian(mu=", round(object@prior@mu,3), ", sigma=", round(object@prior@sigma,3),"), beta=",round(object@beta,3), ", gamma=",round(object@gamma,3)),")\n")
  })
#' @rdname Player
#' @export
performance <- function(a) 0
setGeneric("performance")
#' @rdname Gaussian
#' @export
setMethod("performance", "Player", 
  function(a){
    return(Gaussian(a@prior@mu,sqrt(a@prior@sigma^2 + a@beta^2)))
  })

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

#' @title Game
#'
#' @description Game class
#'
#' @param teams A list of \code{Player} objects. Each position represents a team,
#' so it must contain a vector of \code{Player} objects. 
#' @param result A vector of numbers, with the score obtained by each team, or
#' an empty vector. The default value is an empty vector. In this case, the
#' outcome is defined by the order in which the \code{teams} list was initialized:
#' the teams appearing firstly in the list defeat those appearing later (no ties). If 
#' the list is not empty, it must have the same length as the \code{teams} list. In 
#' this last case, the team with the highest score is the winner, and the teams with 
#' the same score are tied.
#' @param p_draw A number, the probability of a draw. The default value is 
#' \code{P_DRAW = 0}. A rule of thumb states that the probability of a draw must be 
#' initialized with the observed frequency of draws. If in doubt, it is a candidate 
#' parameter to be optimized or integrated by the sum rule. It is used to compute 
#' the prior probability of the observed result, so its value may affect an 
#' eventual model selection task.
#' @param g A game object
#'
#' @return Game object
#'
#' @examples
#' a1 = Player(Gaussian(mu=0, sigma=6), beta=1, gamma=0.03)
#' a2 = Player(); a3 = Player(); a4 = Player()
#' team_a = c(a1, a2)
#' team_b = c(a3, a4)
#' teams = list(team_a, team_b)
#' 
#' g = Game(teams) 
#' post = posteriors(g)
#' lhs = g@likelihoods
#' post[[1]][[1]] == lhs[[1]][[1]]*a1@prior 
#' ev = g@evidence
#' ev == 0.5
#'
#' ta = c(a1)
#' tb = c(a2, a3)
#' tc = c(a4)
#' teams_3 = list(ta, tb, tc)
#' result = c(1, 0, 0)
#' g3 = Game(teams_3, result, p_draw=0.25)
#'
#' @name Game
#' @export
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
partial_evidence <- function(d, margin, tie, e){
    mu = d[[e]]@prior@mu; sigma = d[[e]]@prior@sigma
    return( if (tie[e]) (cdf(margin[e],mu,sigma)-cdf(-margin[e],mu,sigma)) else 1-cdf(margin[e],mu,sigma) )
}

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
    return(list("o"=o, "t"=t, "d"=d, "tie"=tie, "margin"=margin, "evidence" = evidence))
}
#graphical_model = cmpfun(graphical_model )
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
#likelihood_analitico = cmpfun(likelihood_analitico )
likelihood_teams <- function(teams,result,p_draw,gr){
    o= gr$o; d = gr$d; t = gr$t; margin = gr$margin; tie = gr$tie 
    evidence = 1
    step = c(Inf,Inf); i = 0
    while (gr_tuple(step,1e-3) & i < 10){
      for (e in seq(length(d)-1)){
        d[[e]]@prior = posterior_win(t[[e]]) - posterior_lose(t[[e+1]])
        if(i == 0){evidence = evidence * partial_evidence(d, margin, tie, e)}
        d[[e]]@likelihood = approx(d[[e]]@prior,margin[[e]],tie[[e]])/d[[e]]@prior
        likelihood_lose = posterior_win(t[[e]]) - d[[e]]@likelihood
        step = max_tuple(step,delta(t[[e+1]]@likelihood_lose,likelihood_lose))
        t[[e+1]]@likelihood_lose = likelihood_lose
      }
      for (e in seq(length(d),2,-1)){
        d[[e]]@prior = posterior_win(t[[e]]) - posterior_lose(t[[e+1]])
        if((i == 0) & (e == length(d))){evidence = evidence * partial_evidence(d, margin, tie, e)}
        d[[e]]@likelihood = approx(d[[e]]@prior,margin[[e]],tie[[e]])/d[[e]]@prior
        likelihood_win = posterior_lose(t[[e+1]]) + d[[e]]@likelihood
        step = max_tuple(step,delta(t[[e]]@likelihood_win,likelihood_win))
        t[[e]]@likelihood_win = likelihood_win
      }
      i = i + 1
    }
    if (length(d)==1){
      evidence = partial_evidence(gr$d, gr$margin, gr$tie, 1)
      d[[1]]@prior = posterior_win(t[[1]]) - posterior_lose(t[[2]])
      d[[1]]@likelihood = approx(d[[1]]@prior,margin[[1]],tie[[1]])/d[[1]]@prior
    }
    t[[1]]@likelihood_win = posterior_lose(t[[2]]) + d[[1]]@likelihood
    t[[length(t)]]@likelihood_lose = posterior_win(t[[length(t)-1]]) - d[[length(d)]]@likelihood
    res = vector('list', length(t))
    for (e in seq(length(t))){ res[[e]] = likelihood(t[[o[e]]])}
    return(list("messages"=res, "evidence"= evidence))
}
#likelihood_teams = cmpfun(likelihood_teams )
compute_likelihoods <-  function(teams,result,p_draw){
    gr = graphical_model(teams,result,p_draw)
    if (length(teams)>2){
      lhoods = vector('list', length(teams))
      lht = likelihood_teams(teams,result,p_draw,gr)
      m_t_ft = lht$messages
      for (e in seq(length(teams))){#e=2
        lhoods[[e]]  = vector('list', length(teams[[e]]))
        for (i in seq(length(teams[[e]]))){#i=1
          team_perf = N00
          for (a in teams[[e]]){ team_perf  = team_perf + performance(a)}
          ex = exclude(team_perf,teams[[e]][[i]]@prior)
          lhoods[[e]][[i]] = m_t_ft[[e]] - ex
        }
      }
      return(list("likelihoods"=lhoods, "evidence"=lht$evidence))
    }else{
      evidence = partial_evidence(gr$d, gr$margin, gr$tie, 1)
      return(list("likelihoods"=likelihood_analitico(teams,result,p_draw,gr), "evidence"=evidence))
    }
}
#compute_likelihoods = cmpfun(compute_likelihoods)
#' @rdname Game
#' @export
posteriors <- function(g) 0
setGeneric("posteriors")
#' @rdname Game
#' @export
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
          skills[[item$name]]$likelihood <<- skills[[item$name]]$likelihood/item$likelihood * g@likelihoods[[t]][[i]]
          item$likelihood = g@likelihoods[[t]][[i]]
          i = i + 1
        }
        t = t + 1
      }
      events[[e]]$evidence <<- g@evidence
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
      skills[[a]]$backward <<- agents[[a]]$message
    }
    return(iteration())
  },
  new_forward_info = function(){
    for (a in names(skills)){
      skills[[a]]$forward <<- agents[[a]]$receive(skills[[a]]$elapsed)
    }
    return(iteration())
  }
)

#' @title History
#' 
#' @description History class
#' @import hash
#'
#' @param composition A list of list of player's names (id). Each position of the 
#' list is a list that represents the teams of a game, so the latter must contain
#' vectors of names representing the composition of each team in that game.
#' @param results A list of numeric vectors, representing the outcome of each 
#' game. It must have the same 
#' length as the \code{composition} list or be an empty list. The default value is an empty list.
#' When the list is empty, the outcomes of the games are inferred by the order of
#' the teams in the \code{composition} list: the teams appearing firstly in the list 
#' defeat those appearing later (no ties).
#' When the list is not empty, each vector of the list must represents the score of each team in the game. The team with the highest score is 
#' the winner, and the teams with the same score are tied.
#' @param times A numeric vector, the timestamp of each game. It must have the 
#' same length as the \code{composition} list or be an empty list. The default value 
#' is an empty list. When the list is empty, all players' games are 
#' separated by a single timestamp, so a single dynamic uncertainty \code{gamma} will be added between games.
#' When the list is not empty, the amount of dynamic uncertainty will depend on the difference (measured in timestamps) that each player has between games.
#' In addition, the order of the timestamps determines the reading order of the \code{composition} and \code{results} lists.
#' If a player appears more than once in the same timestamp, no dynamic uncertainty will be added between those games.
#' @param priors A hash object, a dictionary of \code{Player} objects indexed by the 
#' players' name (id). Used to create players with special values. The default
#' value is an empty hash. In this case, one \code{Player} object for each unique 
#' name in the \code{composition} list is automatically initialized using the values 
#' of the parameters \code{mu}, \code{sigma}, \code{beta}, and \code{gamma}.
#' The names that appear in the hash are the only ones that will be initialized 
#' with special values.
#' @param mu A number, the prior mean. The deafult value is: \code{MU = 0}.
#' @param sigma A number, the prior standar deviation. The deafult value is: 
#' \code{SIGMA = 6}.
#' @param beta A number, the standard deviation of the performance. The default 
#' value is: \code{BETA = 1}. The parameter \code{beta} acts as the scale of the 
#' estimates. A real difference of one \code{beta} between two skills is equivalent 
#' to 76\% probability of winning.
#' @param gamma A number, the amount of uncertainty (standar deviation) added to 
#' the estimates between events. The default value is: \code{GAMMA = 0.03}.
#' @param p_draw A number, the probability of a draw. The default value is 
#' \code{P_DRAW = 0}. A rule of thumb states that the probability of a draw must be 
#' initialized with the observed frequency of draws. If in doubt, it is a candidate 
#' parameter to be optimized or integrated by the sum rule. It is used to compute 
#' the prior probability of the observed result, so its value may affect an 
#' eventual model selection task.
#' @param epsilon A number, the convergence threshold. Used to stop the convergence procedure. The default value is \code{EPSILON = 1e-6}.
#' @param iterations A number, the maximum number of iterations for convergence. Used to stop the convergence procedure. The default value is \code{ITERATIONS = 30}.
#'
#' @return History object
#' 
#' @field size A number, the amount of games.
#' @field batches A vector of \code{Batch} objects. Where the games that occur at the same timestamp live.
#' @field agents A hash, a dictionary indexed by the players' name (id).
#' @field time A boolean, indicating whether the history was initialized with timestamps or not.
#' @field mu A number, the default prior mean in this particular \code{History} object
#' @field sigma A number, the default prior standard deviation in this particular \code{History} object
#' @field beta A number, the default standar deviation of the performance in this particular \code{History} object
#' @field gamma A number, the default dynamic uncertainty in this particular \code{History} object
#' @field p_draw A number, the probability of a draw in this particular \code{History} object
#' @field h_epsilon A number, the convergence threshold in this particular \code{History} object
#' @field h_iterations A number, the maximum number of iterations for convergence in this particular \code{History} object
#' 
#' @examples 
#' c1 = list(c("a"),c("b"))
#' c2 = list(c("b"),c("c"))
#' c3 = list(c("c"),c("a"))
#' composition = list(c1,c2,c3)
#' h = History(composition, gamma=0.0)
#'
#' trueskill_learning_curves = h$learning_curves()
#' ts_a = trueskill_learning_curves[["a"]]
#' ts_a[[1]]$N; ts_a[[2]]$N
#' ts_a[[1]]$t; ts_a[[2]]$t
#' h$convergence()
#' trueskillThrougTime_learning_curves = h$learning_curves()
#' ttt_a = trueskillThrougTime_learning_curves[["a"]]
#' ttt_a[[1]]$N; ttt_a[[2]]$N
#' ttt_a[[1]]$t; ttt_a[[2]]$t
#' 
#' \dontrun{
#' # Synthetic example
#' library(hash)
#' N = 100
#' skill <- function(experience, middle, maximum, slope){
#' return(maximum/(1+exp(slope*(-experience+middle)))) }
#' target = skill(seq(N), N/2, 2, 0.075)
#' opponents = rnorm(N,target,0.5)
#' composition = list(); results = list(); times = c(); priors = hash()
#' for(i in seq(N)){composition[[i]] = list(c("a"), c(toString(i)))}
#' for(i in
#' seq(N)){results[[i]]=if(rnorm(1,target[i])>rnorm(1,opponents[i])){c(1,0)}else{c(0,1)}}
#' for(i in seq(N)){times = c(times,i)}
#' for(i in seq(N)){priors[[toString(i)]] = Player(Gaussian(opponents[i],0.2))}
#' h = History(composition, results, times, priors, gamma=0.1)
#' h$convergence(); lc_a = h$learning_curves()$a; mu = c()
#' for(tp in lc_a){mu = c(mu,tp[[2]]@mu)}
#' plot(target)
#' lines(mu)
#'
#' # Plotting learning curves
#'
#' # First solve your own example. Here is a dummy one.
#' agents <- c("a", "b", "c", "d", "e")
#' composition <- list()
#' for (i in 1:500) {
#'  who = sample(agents, 2)
#'  composition[[i]] <- list(list(who[1]), list(who[2]))
#' }
#' h <- History(composition = composition, gamma = 0.03, sigma = 1.0)
#' h$convergence(iterations=6)
#'
#' # Then plot some learning curves
#' lc <- h$learning_curves()
#' colors <- c(rgb(0.2,0.2,0.8), rgb(0.2,0.8,0.2), rgb(0.8,0.2,0.2))
#' colors_alpha <- c(rgb(0.2,0.2,0.8,0.2), rgb(0.2,0.8,0.2,0.2), rgb(0.8,0.2,0.2,0.2))
#' plot(0,0, xlim = c(0, 500), ylim = c(-1, 1), xlab = "t", ylab = "skill", type = "n")
#' for (i in 1:3) {
#'   agent <- agents[i]
#'   t <- c(); mu <- c(); sigma <- c()
#'   for(x in lc[[agent]]){
#'     t <- c(t, x$t )
#'     mu <- c(mu, x$N@mu)
#'     sigma <- c(sigma, x$N@sigma)
#'   }
#'   lines(t, mu, col = colors[i], lwd = 2, type = "l")
#'   polygon(c(t, rev(t)), c(mu + sigma, rev(mu - sigma)), col = colors_alpha[i], border = NA)
#' }
#' legend("topright", legend = agents[1:3], col = colors, lwd = 2)
#'
#' }
#' @export History
#' @exportClass History
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
  initialize = function(composition, results=list(), times=c(), priors=hash(),  mu=MU, sigma=SIGMA, beta=BETA, gamma=GAMMA, p_draw=P_DRAW, epsilon=EPSILON,  iterations=ITERATIONS){
    " "
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
        agents[[a]]$last_time <<- if (!time) Inf else t
        agents[[a]]$message <<- b$forward_prior_out(a)
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
        agents[[a]]$message <<- batches[[j+1]]$backward_prior_out(a)
      }
      old = batches[[j]]$posteriors()
      batches[[j]]$new_backward_info()
      step = max_tuple(step,list_diff(old, batches[[j]]$posteriors()))
    }
    #clean(agents)
    for (a in names(agents)){agents[[a]]$message <<- Ninf}
    for (j in seq(2,length(batches))){
      for (a in names(batches[[j-1]]$skills)){
        agents[[a]]$message <<- batches[[j-1]]$forward_prior_out(a)
      }
      old = batches[[j]]$posteriors()
      batches[[j]]$new_forward_info()
      step = max_tuple(step,list_diff(old, batches[[j]]$posteriors()))
    }
    
    } #end ELSE
    return(step)
  },
  convergence = function(epsilon = NA, iterations = NA, verbose=TRUE){
    " "
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
    " "
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
  },
  log_evidence = function(){
    " "
    res = 0
    for(b in batches){
        for(e in b$events){
            res = res + log(e$evidence)
        }
    }
    return(res)
  }
)

#' Print list of Gaussian using the python and julia syntax
#' 
#' @param lc.a List of Gaussians 
#' 
#' @return No return value, print lists of Gaussian using the python and julia syntax
#' 
#' @export
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
