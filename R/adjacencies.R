
#' @export
#' @author Barry Rowlingson, Alison Hale
#' @name voronoi_adjacency
#' @title Find neighbours using Voronoi tessellation.
#' @description Uses Voronoi tessellation to find nearest neighbours and their respective distances.  The neighbours are found using R package deldir.
#' @param data (data frame) with id and x, y columns
#' @param formula a formula of the form id~x+y
#' @param scale scale (divide) coordinates by this factor
#' @param PLOT (boolean) plot neighbour links when TRUE.
#' @return list with components
#'   Adjacencies - adjacency matrix (TRUE means there is a neighbour).
#'   Distances - adjacency matrix with distances between neighbours.
#'   NumNeighbours - total number of neighbours.
#' @importFrom deldir deldir
#' @examples
#' 
#' data=data.frame(i=1:6, x=c(1,2,1,3,4,8),y=c(1,1,5,5,4,8))
#' a=voronoi_adjacency(data,i~x+y)
#' 

voronoi_adjacency = function(data, formula, scale=1, PLOT=FALSE){

    data = makeixy(data, formula, scale)
    
    P=dim(data)[1];  # number of rows

    dd = deldir::deldir(data$x,data$y,suppressMsge=TRUE,plotit=PLOT);  # find adjacencies

    ## create adjacency matrix
    A=matrix(FALSE,P,P);
    A[as.matrix(dd$delsgs[,c("ind1","ind2")])] = TRUE;
    A[as.matrix(dd$delsgs[,c("ind2","ind1")])] = TRUE;

    ## create distance matrix
    D=matrix(NA,P,P);
    D[as.matrix(dd$delsgs[,c("ind1","ind2")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
    D[as.matrix(dd$delsgs[,c("ind2","ind1")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
    
    ## create data frame of results
    N=matrix(colSums(A),P,1); # number of adjacencies for each xy$id
  
    return(list(Adjacencies=A, Distances=D, NumNeighbours=N, ids=data[,"id"], coords=data[,c("x","y")]));
}

makeixy <- function(data, formula, scale){
    m = model.frame(formula, data=data)
    if(ncol(m)!=3){
        stop("incorrect adjacency formula: id~x+y needed")
    }
    names(m)=c("id","x","y")
    m[,2]=m[,2]/scale
    m[,3]=m[,3]/scale
    m
}
    
thresh_adjacency <- function(data, formula, thresh, scale=1){
    data = makeixy(data, formula, scale)
    P = nrow(data)
    D = as.matrix(dist(data[,2:3]))
    D[D>thresh] <- NA
    diag(D) <- NA
    A = !is.na(D)
    N=matrix(colSums(A),P,1)
    return(list(Adjacencies = A, Distances = D, NumNeighbours = N, ids=data[,"id"], coords=data[,c("x","y")]))
}

knn_adjacency <- function(data, formula, k, scale=1){
    data = makeixy(data, formula, scale)
    P = nrow(data)
    nn = as.matrix(reshape2::melt(FNN::knn.index(data[,c("x","y")],k))[,c("Var1","value")])
    A = matrix(FALSE, P, P)
    A[nn]=TRUE
    A[nn[,2:1]] = TRUE
    D = as.matrix(dist(data[,2:3]))
    D[!nn] <- NA
    N=matrix(colSums(A),P,1)
    return(list(Adjacencies = A, Distances = D, NumNeighbours = N, ids=data[,"id"], coords=data[,c("x","y")]))
}

plot_adjacency <- function(adja){
    plot(adja$coords, asp=1)
    dimnames(adja$Adjacencies)=NULL

    segs = reshape2::melt(adja$Adjacencies)
    segs = segs[segs$value,1:2]
    segs = segs[segs$Var1<segs$Var2,]
    x0 = adja$coords[segs[,1],1]
    x1 = adja$coords[segs[,2],1]
    y0 = adja$coords[segs[,1],2]
    y1 = adja$coords[segs[,2],2]
    segments(x0,y0,x1,y1)
}
