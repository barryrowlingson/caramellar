##' Compute exceedences
##'
##' given an array of samples and a first date, compute confidence
##' intervals and exceedences probabilities of given levels and values
##' @title Compute exceedences
##' @param carmcout carmc output
##' @param t1 first date
##' @param clevels confidence interval levels
##' @param elevels exceedence thresholds
##' @return a list of mu [nsites,ntimes] - sample mean; exceedence$e (exceedence levels):
##' exceedence$p[[nlevels]][nsites, ntimes] (exceedence probabilities),
##' cis[[nsites]][ntimes x ncilevels, [min,max,level,date]], dates
##' @author Barry Rowlingson
##' @export
exceeds <- function(carmcout, t1=attr(asamples,"start_date"), clevels=c(.90,.95), elevels=c(2,4,8)){

    asamples = carmcout$S
    ids = attr(asamples,"ids")
    # asamples = mc2arry(samples, nsites, ntimes)

    ## 1. Samples
    ## 2. Times
    ## 3. Sites
    ## asamples[i,j,k] is sample i at time j for site k
    ## will mostly be doing aggregates over c(3,2)

    dates = seq(as.Date(t1), by=1, len=dim(asamples)[3])

    nsites = dim(asamples)[2]
    
    mu = apply(asamples, c(2,3), mean)
    # mu[space, time]
    
    ep <- function(m, e){
        apply(m, 2, function(v){sum(v>e)/length(v)})
    }

    ## apply over sites
    epe <- function(e){
        apply(asamples, 3, function(m){ep(m,e)})
    }

    exa <- function(e){
        apply(asamples, c(3,2), function(v){
            sum(v>e)/length(v)
        }
          )
    }
        
    
    exceeds = list(
        p = lapply(elevels, exa),
        e = elevels
        )

    cis = lapply(1:nsites, function(i){
        do.call(rbind,
                lapply(clevels,function(level){
                    g = (1-level)/2
                    edges = c(g, 1-g) 
                    cs = apply(asamples[,i,], 2, quantile, edges)
                    ms = data.frame(t(cs))
                    names(ms)=c("min","max")
                    ms$level=level
                    ms$date=dates
                    ms
                    
                })
                )
    })
        
    L = structure(list(exceedence = exceeds, mu=mu, cis=cis, dates=dates),
              class="excedence")
    attr(L,"ids") <- ids
    L
}
