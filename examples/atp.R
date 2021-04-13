source("../TrueSkill.R")

data = read.csv("input/history.csv", header=T)

get_composition = function(x){
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

composition =  apply(data, 1, get_composition ) 
times = as.numeric(as.Date(data[,"time_start"], format = "%Y-%m-%d"))
#as.Date(times[length(times)], origin=as.Date("1970-01-01", format = "%Y-%m-%d"))

start.time <- Sys.time()
h = History(composition = composition, times = times, sigma = 1.6, gamma = 0.036)
h$convergence(epsilon=0.01, iterations=10)
end.time <- Sys.time()
print(end.time - start.time)



