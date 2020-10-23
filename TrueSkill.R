library(microbenchmark)

BETA = 1.0
MU = 0.0
SIGMA = BETA * 6
GAMMA = BETA * 0.05
P_DRAW = 0.0
EPSILON = 1e-6
ITERATIONS = 10
PI = 1/(SIGMA^2)
TAU = MU*PI

Environment <- setRefClass(
  "Environment", 
  fields = list(mu = "numeric"
               ,sigma = "numeric"
               ,beta = "numeric"
               ,gamma = "numeric"
               ,p_draw = "numeric"
               ,epsilon = "numeric"
               ,iterations = "numeric")
)
Environment$methods( 
  initialize = function(mu=MU, sigma=SIGMA, beta=BETA, gamma=GAMMA, p_draw=P_DRAW, epsilon=EPSILON, iterations=ITERATIONS) {
    mu <<- mu; sigma <<- sigma; beta <<- beta; gamma <<-gamma; p_draw <<- p_draw; epsilon <<- epsilon; iterations <<- iterations
    },
  show = function() { cat(paste0("Environment(mu=", round(mu,3), ", sigma=", round(sigma,3),", beta=",round(beta,3),", gamma=",round(gamma,3),", p_draw=", round(p_draw),", epsilon=",epsilon,", iterations=",ceiling(iterations),")"))}
)

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

cdf = function(x, mu, sigma){
  z = -(x - mu) / (sigma * sqrt(2))
  return(0.5*erfc(z) )
}
pdf = function(x, mu, sigma){
  normalizer = (sqrt(2*pi) * sigma)^-1
  functional = exp( -((x - mu)^2) / (2*sigma ^2) ) 
  return(normalizer * functional)
}
ppf = function(p, mu, sigma){
  return( mu - sigma * sqrt(2)  * erfcinv(2 * p))
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
  return(Gaussian(mu_trunc, sigma_trunc))
}
compute_margin = function(p_draw, sd){
  return(abs(ppf(0.5-p_draw/2, 0.0, sd )))
}
max_tuple = function(t1, t2){
  return(c(max(t1[1],t2[1]), max(t1[2],t2[2])))
}
gr_tuple = function(tup, threshold){
    return( (tup[1] > threshold) | (tup[2] > threshold) )
}
sortperm = function(xs){
  return(order(xs))
}
Gaussian <- setRefClass("Gaussian", 
                fields = list(mu = "numeric"
                             ,sigma = "numeric"))
Gaussian$methods( 
  initialize = function(mu=MU,sigma=SIGMA) {
    if(sigma>=0.0){
      mu <<- mu; sigma <<- sigma
    }else{
      stop("Require: (sigma >= 0.0)")
    }
  },
  show = function() { cat(paste0("Gaussian(mu=", round(mu,3), ", sigma=", round(sigma,3),")"))},
  delta = function(M){
    return( c(abs(mu - M$mu) , abs(sigma - M$sigma)) )
  },
  pi = function(){
    if (sigma > 0.0){return(sigma^-2)
    }else{return(Inf)}
  },
  tau = function(){
    if (sigma > 0.0){return(mu*sigma^-2)
    }else{return(Inf)}
  },
  forget = function(gamma, t){
    return(Gaussian(mu,sqrt(sigma^2 + t*gamma^2)))
  },
  exclude = function(M){
    return(Gaussian(mu-M$mu, sqrt(sigma^2 - M$sigma^2)))
  },
  isapprox = function(M,tol=1e-4){
    return(abs(mu-M$mu)<tol & abs(sigma-M$sigma)<tol)
  }
)
setMethod("+", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    mu = e1$mu + e2$mu
    sigma = sqrt((e1$sigma^2) + (e2$sigma^2) )
    return(Gaussian(mu,sigma))
})
setMethod("-", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    mu = e1$mu - e2$mu
    sigma = sqrt((e1$sigma^2) + (e2$sigma^2) )
    return(Gaussian(mu,sigma))
})
setMethod("*", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    tau_ = e1$tau() + e2$tau(); pi_ =  e1$pi() + e2$pi()
    mu.sigma = mu_sigma(tau_, pi_)
    return(Gaussian(mu.sigma[1],mu.sigma[2]))
})
setMethod("/", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    tau_ = e1$tau() - e2$tau(); pi_ =  e1$pi() - e2$pi()
    mu.sigma = mu_sigma(tau_, pi_)
    return(Gaussian(mu.sigma[1],mu.sigma[2]))
})
setMethod("==", c("Gaussian", "Gaussian"),
  function(e1, e2) {
    mu = e1$mu - e2$mu < 1e-3
    sigma = if (e2$sigma == Inf | e1$sigma == Inf) (e1$sigma==e2$sigma) else  (e1$sigma - e2$sigma < 1e-3)
    return(mu & sigma)
})


N01 = Gaussian(0,1)
N00 = Gaussian(0,0)
Ninf = Gaussian(0,Inf)
Nms = Gaussian(MU, SIGMA)

list_diff = function(old, new){
  step = c(0,0)
  for (a in names(old)){#a="0"
    step = max_tuple(step, old[[a]]$delta(new[[a]]))
  }
  return(step)
}


Rating <- setRefClass("Rating", 
  fields = list(N = "Gaussian"
               ,beta = "numeric"
               ,gamma = "numeric"
               ,draw = "Gaussian")
)
Rating$methods( 
  initialize = function(mu=MU, sigma=SIGMA, beta=BETA, gamma=GAMMA){
    N <<- Gaussian(mu, sigma); beta <<- beta; gamma <<- gamma; draw <<- Ninf
  },
  show = function() { cat(paste0("Rating(mu=", round(N$mu,3), ", sigma=", round(N$sigma,3),")"))},
  performance = function(){
    return(Gaussian(N$mu, sqrt(N$sigma^2 + beta^2)))
  }
)

team_messages <- setRefClass("team_messages",
  fields = list(
    prior = "Gaussian",
    likelihood_lose = "Gaussian",
    likelihood_win = "Gaussian",
    likelihood_draw = "Gaussian")
)
team_messages$methods(
  initialize = function(prior=Ninf, likelihood_lose=Ninf, likelihood_win=Ninf, likelihood_draw=Ninf){
    prior <<- prior;
    likelihood_lose <<- likelihood_lose 
    likelihood_win <<- likelihood_win
    likelihood_draw <<- likelihood_draw
  },
  p = function(){
    return(prior*likelihood_lose*likelihood_win*likelihood_draw)
  },
  posterior = function(){
    return(prior*likelihood_lose*likelihood_draw)
  },
  posterior_win = function(){
    return(prior*likelihood_lose*likelihood_draw)
  },
  posterior_lose = function(){
    return(prior*likelihood_win*likelihood_draw)
  },
  likelihood = function(){
    return(likelihood_win*likelihood_lose*likelihood_draw)
  }
)

draw_messages <- setRefClass("draw_messages",
  fields = list(
    prior = "Gaussian",
    prior_team = "Gaussian",
    likelihood_lose = "Gaussian",
    likelihood_win = "Gaussian")
)
draw_messages$methods(
  initialize = function(prior=Ninf, prior_team=Ninf, likelihood_lose=Ninf, likelihood_win=Ninf){
    prior <<- prior;
    prior_team <<- prior_team
    likelihood_lose <<- likelihood_lose
    likelihood_win <<- likelihood_win
  },
  p = function(){
    return(prior_team*likelihood_lose*likelihood_win)
  },
  posterior_win = function(){
    return(prior_team*likelihood_lose)
  },
  posterior_lose = function(){
    return(prior_team*likelihood_win)
  },
  likelihood = function(){
    return(likelihood_win*likelihood_lose)
  }
)

diff_messages <- setRefClass("diff_messages ",
  fields = list(
    prior = "Gaussian",
    likelihood = "Gaussian")
)
diff_messages$methods(
  initialize = function(prior=Ninf, likelihood=Ninf){
    prior <<- prior
    likelihood <<- likelihood
  },
  p = function(){
    return(prior*likelihood)
  }
)

Game <- setRefClass("Game",
  fields = list(
    teams = "list",
    result = "vector",
    p_draw = "numeric",
    likelihoods = "list",
    evidence = "numeric")
)
Game$methods(
  initialize = function(teams, result, p_draw=0.0){
    if (length(teams) != length(result)) stop("length(teams) != length(result)")
    if ((0.0 > p_draw) | (1.0 <= p_draw)) stop("0.0 <= p_draw < 1.0")
        
    teams <<- teams
    result <<- result
    p_draw <<- p_draw
    likelihoods <<- list()
    evidence <<- 0.0
    compute_likelihoods()
  },
  size = function(){
    res = c()
    for (team in teams){ res = c(res,length(team)) }
    return(res)
  },
  performance = function(i){
    res = Gaussian(0,0)
    for (r in teams[[i]]){ res = res + r$performance()}
    return(res)
  }
  ,
  graphical_model = function(){
    o = sortperm(result)
    t = c(); d = c(); tie = c(); margin = c()
    for (e in seq(length(teams))){
      t = c(t,team_messages(performance(o[e]))) }
    for (e in seq(length(teams)-1)){ 
      d = c(d,diff_messages(t[[e]]$prior - t[[e+1]]$prior)) }
    for (e in seq(length(d))){
     tie = c(tie,result[o[e]]==result[o[e+1]])}
    for (e in seq(length(d))){
      if (p_draw == 0.0){ 
        margin = c(margin, 0.0)}
      else{ 
        betas = 0
        for (a in teams[[o[e]]]){betas = betas + a$beta^2}
        for (a in teams[[o[e+1]]]){betas = betas + a$beta^2}
        margin = c(margin, compute_margin(p_draw, sqrt(betas) )) 
      }
    }
    evidence <<- 1
    for (e in seq(length(d))){
        mu = d[[e]]$prior$mu; sigma = d[[e]]$prior$sigma
        evidence <<- evidence * (if (tie[e]) (cdf(margin[e],mu,sigma)-cdf(-margin[e],mu,sigma)) else 1-cdf(margin[e],mu,sigma))
    }
    return(list("o"=o, "t"=t, "d"=d, "tie"=tie, "margin"=margin))
  },
  likelihood_analitico = function(){
    gr = graphical_model()
    d = gr$d[[1]]$prior
    mu_sigma_trunc =  trunc(d$mu, d$sigma, gr$margin[1], gr$tie[1])
    mu_trunc = mu_sigma_trunc$mu; sigma_trunc = mu_sigma_trunc$sigma
    if (d$sigma==sigma_trunc){
      delta_div = 0.0
      theta_div_pow2 = Inf
    }else{
      delta_div = (d$sigma^2*mu_trunc - sigma_trunc^2*d$mu)/(d$sigma^2-sigma_trunc^2)
      theta_div_pow2 = (sigma_trunc^2*d$sigma^2)/(d$sigma^2 - sigma_trunc^2)
    }
    res = list()
    for (i in seq(length(teams))){#i=1
      team = c()
      for (j in seq(length(teams[[i]]))){#j=1
        mu = if (d$sigma == sigma_trunc) 0.0 else teams[[i]][[j]]$N$mu + ( delta_div - d$mu)*(-1)^(gr$o[i]==2)
        sigma_analitico = sqrt(theta_div_pow2 + d$sigma^2 - teams[[i]][[j]]$N$sigma^2)
        team = c(team,Gaussian(mu,sigma_analitico))
      res[[i]] = team
      }
    }
    return(res)
  },
  likelihood_teams = function(){
    gr = graphical_model()
    o= gr$o; d = gr$d; t = gr$t; margin = gr$margin; tie = gr$tie 
    step = c(Inf,Inf); i = 0
    while (gr_tuple(step,1e-6) & i < 10){
      for (e in seq(length(d)-1)){
        d[[e]]$prior = t[[e]]$posterior_win() - t[[e+1]]$posterior_lose()
        d[[e]]$likelihood = trunc(d[[e]]$prior$mu,d[[e]]$prior$sigma,margin[[e]],tie[[e]])/d[[e]]$prior
        likelihood_lose = t[[e]]$posterior_win() - d[[e]]$likelihood
        step = max_tuple(step,t[[e+1]]$likelihood_lose$delta(likelihood_lose))
        t[[e+1]]$likelihood_lose = likelihood_lose
      }
      for (e in seq(length(d),2,-1)){
        d[[e]]$prior = t[[e]]$posterior_win() - t[[e+1]]$posterior_lose()
        d[[e]]$likelihood = trunc(d[[e]]$prior$mu,d[[e]]$prior$sigma,margin[[e]],tie[[e]])/d[[e]]$prior
        likelihood_win = t[[e+1]]$posterior_lose() + d[[e]]$likelihood
        step = max_tuple(step,t[[e]]$likelihood_win$delta(likelihood_win))
        t[[e]]$likelihood_win = likelihood_win
      }
      i = i + 1
    }
    if (length(d)==1){
      d[[1]]$prior = t[[1]]$posterior_win() - t[[2]]$posterior_lose()
      d[[1]]$likelihood = trunc(d[[1]]$prior$mu,d[[1]]$prior$sigma,margin[[1]],tie[[1]])/d[[1]]$prior
    }
    t[[1]]$likelihood_win = t[[2]]$posterior_lose() + d[[1]]$likelihood
    t[[length(t)]]$likelihood_lose = t[[length(t)-1]]$posterior_win() - d[[length(d)]]$likelihood
    res = c()
    for (e in seq(length(t))){ res = c(res,t[[o[e]]]$likelihood())}
    return(res)
  },
  compute_likelihoods = function(){
    if (length(teams)>2){
      m_t_ft = likelihood_teams()
      for (e in seq(length(teams))){#e=1
        res = c()
        for (i in seq(length(teams[[e]]))){#i=1
          res = c(res, m_t_ft[[e]] - performance(e)$exclude(teams[[e]][[i]]$N))
        }
        likelihoods[[e]] <<- res
      }
    }else{
      likelihoods <<- likelihood_analitico()
    }
  },
  posteriors = function(){
    res = list()
    for (e in seq(length(teams))){
      post = c()
      for (i in seq(length(teams[[e]]))){
        post = c(post, likelihoods[[e]][[i]] * teams[[e]][[i]]$N)   
      }
      if (is.null(names(teams))){
        res[[e]] = post
      }else{
        res[[names(teams)[e]]] = post
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
    return(forward * likelihood * backward)
  },
  posterior_back = function(){
    return(forward * likelihood)
  },
  posterior_for = function(){
    return(likelihood * backward)
  }
)

Agent <- setRefClass("Agent",
  fields = list(
    prior = "Rating",
    message = "Gaussian",
    last_time = "numeric")
)
Agent$methods(
  initialize = function(prior, message, last_time){
    prior <<- prior; message <<- message; last_time <<- last_time
  },
  receive = function(elapsed){
    if (!(message==Ninf)){
      res = message$forget(prior$gamma, elapsed)
    }else{
      res = prior$N
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
  show = function(){
    cat(paste0("Event(",names(),result()))
  },
  names = function(){
    res = list()
    for (t in seq(length(teams))){
      vec = c()
      for (item in teams[t]$items){
        vec = c(vec, item$name)
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
  res = c()
  for (xs in xss){for (x in xs){res = c(res, x)}}
  return(unique(res))
}

initialize_events = function(composition, results){
    events_ = c()
    for (e in seq(length(composition))){
      teams_ = c()
      for (t in seq(length(composition[[e]]))){
        items_ = c()
        for (a in seq(length(composition[[e]][[t]]))){
          items_ = c(items_, Item(composition[[e]][[t]][[a]],Ninf))
        }
        teams_ = c(teams_, Team(items_,results[[e]][[t]]))
      }
      events_ = c(events_, Event(teams_, 0))
    }
    return(events_)
}
initialize_skills = function(composition,agents,time){
    this_agents = list_unique(composition)
    skills_ = list()
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
    skills = "list",
    agents = "list",
    env = "Environment"
    )
)
Batch$methods(
  initialize = function(composition, results ,time, agents, env){
    if (length(composition) != length(results)) stop("length(composition)!= length(results)")
      
    skills <<- initialize_skills(composition, agents, time)
    events <<- initialize_events(composition, results)
    time <<- time
    agents <<- agents
    env <<- env
    iteration()
  },
  show = function(){
    cat(paste0("Batch(time=",time,", events=",length(events),", skills=", length(skills),")\n"))
  },
  posterior = function(a){
    return(skills[[a]]$posterior() )
  },
  posteriors = function(){
    res = list()
    for (a in names(skills)){
      res[[a]] = skills[[a]]$posterior() 
    }
    return(res)
  },
  within_prior = function(item){
    prior = agents[[item$name]]$prior
    ms = skills[[item$name]]$posterior()/item$likelihood
    return(Rating(ms$mu, ms$sigma, prior$beta, prior$gamma))
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
      g = Game(teams, result, env$p_draw)
      t = 1
      for (team in events[[e]]$teams){
        i = 1
        for (item in team$items){
          skills[[item$name]]$likelihood <<- (skills[[item$name]]$likelihood / item$likelihood) * g$likelihoods[[t]][[i]]
          item$likelihood = g$likelihoods[[t]][[i]]
          i = i + 1
        }
        t = t + 1
      }
      events[[e]]$evidence <<- g$evidence
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
    return(skills[[name]]$posterior_for()$forget(agents[[name]]$prior$gamma,skills[[name]]$elapsed))
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


History = setRefClass("History",
  fields = list(
    size = "numeric",
    batches = "vector",
    agents = "list",
    time = "logical",
    env = "Environment")
)
History$methods(
  initialize = function(composition,results,times=c(),priors=list(), env=Environment()){
    if (length(composition) != length(results)){ stop("length(composition) != length(results)")}
    if (length(times) > 0 & (length(composition) != length(times))){ stop("length(times) error")}
    
    this_agents = list_unique(composition)
    agents_ = list()
    for (a in this_agents ){
        agents_[[a]] = Agent(if (a %in% names(priors)) priors[[a]] else Rating(env$mu, env$sigma, env$beta, env$gamma), Ninf, -Inf)
    }
    
     size <<- length(composition)
     agents <<- agents_
     env <<- env
     time <<- length(times) > 0
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
      b = Batch(composition[o[seq(i,j)]], results[o[seq(i,j)]], t, agents, env)
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
  convergence = function(verbose=FALSE){
    step = c(Inf, Inf); i = 1
    while (gr_tuple(step,env$epsilon) & i <= env$iterations){
      if (verbose){cat(paste0("Iteration = ", i))}
      step = iteration()
      i = i + 1
      if (verbose){cat(paste0(" step = (", step[1], ", ", step[2], ")\n"))}
    }
    if (verbose){cat("End\n")}
    return(step)
  }
)









# 
# teams = list(ta = c(Rating()), tb = c(Rating()))
# result = c(0,1)
# g = Game(teams, result)
# g$size()
# microbenchmark(Rating(0.001,1.0), times=5, unit="s")
# microbenchmark(N1/N2, times=5, unit="s")
# 
# 
