#' @export
plot_betas <- function(stout, type="line", split=1, iteration=FALSE){

    melted = reshape2::melt(stout$beta, varnames=c("i","Beta"))
    if(iteration){
        melted$i = stout$iteration(melted$i)
        xl <- "Iteration"
    }else{
        xl <- "Sample"
    }
    if(split>1){
        maxI <- max(melted$i)
        minI <- min(melted$i)
        sVal <- (melted$i-minI)/(maxI-minI)
        melted$part <- paste("Part",1+trunc(split * (sVal) - 0.001))
    }
    if(type=="line"){
        return(
        ggplot2::ggplot(
            melted,
            aes(x=i, y=value,group=Beta)) +
                ggplot2::geom_line() +
                    ggplot2::facet_wrap(~Beta,scales="free_y") +
            xlab(xl)
            )
    }

    if(type=="hist"){
        if(split>1){
            f=~Beta+part
        }else{
            f=~Beta
        }
        return(
            ggplot2::ggplot(melted, aes(x=value)) + geom_histogram() + facet_wrap(f, scales="free")
            )
    }
    if(type=="boxplot"){
        if(split>1){
            g = ggplot(melted,aes(x=part, y=value))
        }else{ 
            g = ggplot(melted,aes(x="", y=value))
        }
            
        return(g+geom_boxplot() + facet_wrap(~Beta, nrow=1, scales="free_y") +
               theme(axis.ticks.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank())
               )
    }
}

#' @export
plot_vars <- function(stout, type="line", split=1, iteration=FALSE){
    if(iteration){
        xl <- "Iteration"
        i = stout$iteration(1:length(stout$delta))
    }else{
        xl <- "Sample"
        i = 1:length(stout$delta)
    }
    d = data.frame(i=i,delta=stout$delta, sigma2=stout$sigma2, phi=stout$phi, rho=stout$rho)
    d = reshape2::melt(d, id.vars="i")

    if(split>1){
        maxI <- max(d$i)
        minI <- min(d$i)
        sVal <- (d$i-minI)/(maxI-minI)
        d$part <- paste("Part",1+trunc(split * (sVal) - 0.001))
    }

    
    if(type=="line"){
        return(ggplot(d, aes(x=i, y=value)) + geom_line() + facet_wrap(~variable,ncol=1,scales="free_y") + xlab(xl))
    }
    if(type=="hist"){
        if(split>1){
            g = ggplot(d, aes(x=value)) + geom_histogram() + facet_grid(part~variable, scales="free")
        }else{
            g = ggplot(d, aes(x=value)) + geom_histogram() + facet_wrap(~variable,scales="free_x")
        }
        return(g)
    }
    if(type=="boxplot"){
        if(split<=1){
            return(ggplot(d,aes(x="", y=value))+geom_boxplot() + facet_wrap(~variable, nrow=1, scales="free_y") +
                   theme(axis.ticks.x=element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank())
                   )
        }
        return(ggplot(d, aes(x=part, y=value)) + geom_boxplot() + facet_wrap(~variable, nrow=1, scales="free_y"))
    }
}

#' @export
plot_exceeds <- function( exceedences, site_id, ev=1, traffic=FALSE, title=""){

    site = location(exceedences, site_id)
    
    d = data.frame(date=exceedences$dates,
        variable="mean",
        rate=exceedences$mu[site,],
        stringsAsFactors=FALSE)

    elevel = exceedences$exceedence$e[ev]
    p = exceedences$exceedence$p[[ev]][,site]
    pCat = cut(p, breaks=c(-1,.8,.9,.99,2),labels=c("Low","Medium","High","VHigh"))
    expdata = data.frame(
        date=exceedences$dates,
        y = 1,
        p = pCat
        )
    
    g = ggplot(data=d, aes(x=date, y=rate)) + geom_line() +
        geom_ribbon(data=exceedences$cis[[site]],aes(x=date,ymax=max,ymin=min,y=min,group=level),alpha=0.2)+ 
        ggtitle(title) + theme(plot.title=element_text(hjust=0))
    if(traffic){
        g = g + geom_point(data=expdata, aes(x=date,y=y,col=p), size=5) +
            scale_colour_manual(values=c("Low"="green","Medium"="yellow","High"="orange","VHigh"="red"),
                                limits=c("Low","Medium","High","VHigh"),
                                name=paste0("P(rate > ",elevel,")"))
    }
    g
}

#' @export
plot.carmc <- function(x,y,...){
    betas <- plot_betas(x,...)
    vars <- plot_vars(x,...)
    grid.arrange(betas,vars)
}

#' @export
plot_S <- function(out, ids, samples=1:nrow(out$S), alpha=0.1, col="#000000", ncol=1, force=FALSE){
    dates = carmcdates(out)
    complexity = length(ids)*length(samples)*length(dates)
    if(complexity > 10*1000*20 & !force){
        stop("Seems like a lot of graphics. Use the force=TRUE parameter if you really want to do this.")
    }
    inds = location(out$S,ids)
    if(any(is.na(inds))){
        stop()
    }
    Sv = melt(out$S[samples,inds,,drop=FALSE])
    names(Sv)=c("Sample","Location","Date","S")
    Sv$Location = ids[Sv$Location]
    Sv$Date = dates[Sv$Date]
    ggplot(Sv, aes(x=Date, y=S, group=Sample)) + geom_line(col=col, alpha=alpha) + facet_wrap(~Location, ncol=ncol)
}

#' @export
acf_vars <- function(out,...){
    d = data.frame(delta=out$delta, sigma2=out$sigma2, phi=out$phi, rho=out$rho)
    qwraps2::qacf(d,...)
}

#' @export
acf_beta <- function(out,...){
    qwraps2::qacf(out$beta,...)
}

#' @export
plot_prop <- function(out,vars=names(out$prop)){
    hsamples = seq(1, out$mc_control$n.sim, by=out$mc_control$hthin)
    d = rbind(
        do.call(rbind,
                lapply(vars, function(var){data.frame(iteration=hsamples, value=out$prop[[var]][[1]], var=var, what="bw")})
                ),
        do.call(rbind,
                lapply(vars, function(var){data.frame(iteration=hsamples, value=out$prop[[var]]$a, var=var, what="acc")})
                )
        )

    ggplot(d, aes(x=iteration, y=value, group=var)) + geom_line() + facet_wrap(~what+var, nrow=2,scales="free_y") 

}
