### code for generating all plots for a single analysis as PNG file in a folder

#' @export
all_plots <- function(out, e, root){
    p = make_image_path(root)
}

ts_plots <- function(out, e, p, width=500, height=250, ...){
    for(id in location_ids(out)){
        path = p(id, "exceeds", "exceed")
        png(path, width=width, height=height)
        print(plot_exceeds(e, id, traffic=TRUE))
        dev.off()              
    }
    
}

#' @export
map_plots <- function(...){

}

#' @export
country_plot <- function(bg){
    if(inherits(bg, "Raster")){
        mytiles = bg
    }else{
        mytiles = mapmisc::openmap(spTransform(bg,CRS("+init=epsg:4326")))
    }
    mapmisc::map.new(mytiles)
    plot(mytiles)
}

#' @export
make_image_path <- function(root){
    force(root)
    p = function(id, dir, name){
        filename = paste0(name,"_",id,".png")
        path = file.path(root, dir, filename)
        dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
        path
    }
    p
}

#' @export
augment_catchments <- function(case_ids, catchments){
    nCases = as.data.frame(table(case_ids))
    names(nCases)=c("id","count")
    catchments@data %<>% dplyr::left_join(nCases)
    catchments$density = catchments$count*(catchments$pc/100)/catchments$area
    catchments
}
