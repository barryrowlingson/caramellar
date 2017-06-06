## compute catchments

##' compute catchment at various percentages
##'
##' given a few percentages, work out the catchment at each level
##' @title multiple catchments
##' @param locations spatial points
##' @param pcs vector of percentages
##' @return spatial polygons with percentage and source id as attribute
##' @author Barry Rowlingson
##' @export
multicatch <- function(locs,pcs, idcolumn=1){


    locs = locs[,idcolumn]

    ## now a bit of cheeky augmentation so all locs have >5 pts
    idcount = as.data.frame(table(locs@data[,1]))
    if(any(idcount$Freq)==0){
        stop("some premises have zero counts")
    }
    smallptcounts = idcount$Var1[idcount$Freq<=5]
    if(length(smallptcounts)>0){
        extras = makeextra(locs, smallptcounts)
        locs = rbind(locs, extras)
    }
    
    ## have to call catchments with a SPDF with one column...
    do.call(rbind.SpatialPolygonsDataFrame,
            lapply(pcs,function(pc){
                cc = catchments(locs,pc);
                cc$pc=pc;
                spChFIDs(cc)=paste0(cc$id,"-",pc)
                cc
            })
            )
}

##' compute catchments
##'
##' use habitat models for catchments
##' @title compute catchments
##' @param locs spatial coordinates
##' @param pc return catchment area at this percentage
##' @return catchment
##' @author Barry Rowlingson
##' @export
catchments <- function(locs, pc){
    adehabitatHR::mcp(locs, percent=pc)
}

makeextra <- function(locs, ids){
    augments = do.call(rbind.SpatialPointsDataFrame,
        lapply(ids, function(id){
            makefakes(locs[locs@data[,1]==id,])
        })
        )
    augments
}

makefakes <- function(locs, include_srcs=FALSE){
    xy = coordinates(locs)
    ## c() here is necessary since replicate returns a matrix
    x = c(replicate(5, jitter(xy[,1])))
    if(include_srcs){
        x = c(xy[,1],x)
    }
    y = c(replicate(5, jitter(xy[,2])))
    if(include_srcs){
        y = c(xy[,2],y)
    }
    newpts = SpatialPointsDataFrame(cbind(x,y),data.frame(Z=rep(locs@data[1,1], length(x))))
    names(newpts)=names(locs)
    newpts
}

    
getxyweights <- function(xy, catchments){
    overs = over(xy,as(catchments,"SpatialPolygons"),returnList=TRUE)[[1]]
    wts = catchments@data[overs,] %>% dplyr::group_by(id) %>% dplyr::filter(pc==min(pc)) %>% dplyr::select(id, density) %>%
        group_by() %>% mutate(p = density/sum(density))
    wts    
}

getxysamples <- function(xy, asamples, tindex, catchments, n){
    wts = getxyweights(xy, catchments)
    if(nrow(wts)==0){
        return(rep(NA,n))
    }
    wtsample(asamples[,,tindex], attr(asamples,"ids"), wts, n)
}

wtsample <- function(samplematrix, idlookup, wts, n){
    pindex = match(wts$id, idlookup)
    npolys = nrow(wts)
    nsamples = nrow(samplematrix)
    samples = samplematrix[,pindex]
    ps = rep(wts$p, rep(nsamples, npolys))
    sample(samples, size=n, p = ps)
}

rasterExceedence <- function(r, thresh, asamples, tindex, catchments, n, progress=FALSE){
    xys = as(r, "SpatialPoints")
    proj4string(xys)=proj4string(catchments)
    if(progress){
        pb = txtProgressBar(min=1, max=length(r), style=3)
    }
    r[] = sapply(seq_along(r), function(i){
        if(progress){setTxtProgressBar(pb,i)}
        sum(getxysamples(xys[i,], asamples, tindex, catchments, n)>thresh)/n
    }, simplify=TRUE)

    r
}

