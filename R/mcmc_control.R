##' mcmc control parameters
##'
##' mcmc run control parameters
##' @title mcmc control parameters
##' @param n.sim number of steps
##' @param burnin number of burn-in steps
##' @param thin sample thinning
##' @param hthin thinning of proposal parameters/acceptance sampling
##' @param tmax maximum run duration, not including burn-in.
##' @return a list of parameters
##' @author Barry Rowlingson
##' @export
mcmc_control = function(n.sim, burnin, thin, hthin, tmax=lubridate::dseconds(Inf)){
    if(!inherits(tmax,"Duration")){
        warning("tmax should be a duration. Attempting to convert...")
        tmax = as.duration(tmax)
    }
    C = list(n.sim=n.sim, burnin=burnin, thin=thin, hthin=thin, tmax=tmax)
    class(C) <- "mcmc_control"
    C
}
