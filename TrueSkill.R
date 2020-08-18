Gaussian <- setRefClass(
  "fooo", 
  fields = list(mu = "numeric"
               ,sigma = "numeric"
               ,tau = "numeric"
               ,pi = "numeric")
)
Gaussian$methods( 
  initialize = function(a=25,b=25/6,inverse=F) {
    if(!inverse){
      mu_ = a; sigma_ = b
      mu <<- mu_; sigma <<- sigma_
      pi <<- sigma_^-2 ; tau <<- mu_ * pi
    }
    if(inverse){
      tau_ = a; pi_ = b
      pi <<- pi_; tau <<- tau_
      mu <<- tau_/pi_ ; sigma <<- sqrt(1/pi_)
    }
  },
  show = function() { cat(paste0("Gaussian(mu=", mu, ", sigma=", sigma,")"))}
)

setMethod("+", c("fooo", "fooo"),
  function(e1, e2) {
    mu_ = e1$mu + e2$mu
    sigma_ = sqrt((e1$sigma^2) + (e2$sigma^2) )
    return(Gaussian(mu_,sigma_))
})

