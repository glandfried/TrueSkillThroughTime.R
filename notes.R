

setClass("Gaussian", representation(mu = "numeric", sigma = "numeric"))
Gaussian <- function(mu=MU,sigma=SIGMA) {
  if(sigma>=0.0){
    return(new("Gaussian",mu=mu,sigma=sigma))
  }else{
    stop("Require: (sigma >= 0.0)")
  }
}
setMethod("+", c("Gaussian", "Gaussian"),
    function(g1, g2) {
        mu = g1@age + g2@age
        sigma = sqrt((g1@sigma^2) + (g2@sigma^2) )
        return(Gaussian(mu,sigma))}
)
Gaussian$methods( 
  show = function() { cat(paste0("Gaussian(mu=", round(mu,3), ", sigma=", round(sigma,3),")"))}
)


Gaussian(0.001,1.0)
microbenchmark(Gaussian(0.001,1.0), times=5, unit="s")

object.size(Gaussian(0.001,1))

