#' @title cluster compactness measures  
#
#' @description Calculate 5 within cluster compactness or homogenerity measures 
#' 
#' @param distance  dissimilarity matrix 
#' @param clusters  existing cluster partitions to determine compactness 
#' @param data  data matrix from which cluster partitions were derived. Note that one of  \code{\link{distance}} 
#' or \code{\link{data}} must be provided
#' @param ensemble  whether to find the optimal compactness measure by majority rule 
#' @param compact.measure if ensemble=FALSE, provide the compactness measure to return. Options are 
#' \itemize{
#' \item DM  cluster diabmeter 
#' \item SW  Cluster Silhouette width
#' \item WD  Within cluster distance 
#' \item gap  Cluster gap
#' } 
#' @param \dots further arguments passed to or from other methods.
#' @return the numeric compactness measure   
#' @export
clusterCompactness <- function(distance=NULL, clusters, Data=NULL, method="euclidean", ensemble = TRUE, compact.measure = "SW", ...) { 

 if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
 if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
 if (class(distance)=="dist") distance <- as.matrix(distance)

mod <- clusterStats(d = distance, clustering = clusters, alt.clustering=NULL, ...)

DM <-  which.min(mod$diameter) ## diameter, small better 
SW <-  which.max(mod$clus.avg.silwidths) # large better 
WD <-  which.min(mod$within.dist) ### mean within distance  - small better  
gap <- which.min(mod$cwidegap)  ### with cluster gap .. small better 
##dia <- which.min(mod$diameter)  ### within cluster diameter ... small better 
SSW <- which.min(mod$within.cluster.ss) ### within cluster variance ... small beter 
tb <- c(DM, SW,WD,gap,SSW) 
if(ensemble) return(getMode(tb))
else {
switch(compact.measure, 
DM ={return(DM)},
SW ={return(SW)},
SSW={return(SSW)},
WD={return(WD)},
gap={return(gap)},
stop("Compact measure choice unknown")
)
}  
}



clusterCompactness.values <- function(distance=NULL, clusters, Data=NULL, method="euclidean", compact.measure = "SW", ...) { 

 if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
 if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
 if (class(distance)=="dist") distance <- as.matrix(distance)

mod <- clusterStats(d = distance, clustering = clusters, alt.clustering=NULL, ...)

switch(compact.measure, 
DM ={return(-1*mod$diameter)},
SW ={return(mod$clus.avg.silwidths)},
SSW={return(-1*mod$within.dist)},
WD={return(-1*mod$cwidegap)},
gap={return(-1*mod$within.cluster.ss)},
stop("Compact measure choice unknown")
)
}



clusterCompactness.ind <- function(distance=NULL, clusters, Data=NULL, method="euclidean", ...) { 

 if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
 if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
 if (class(distance)=="dist") distance <- as.matrix(distance)

mod <- clusterStats(d = distance, clustering = clusters, alt.clustering=NULL, ...)
DM <- rank(-1*mod$diameter, ties.method = "first")## diameter, small better 
SW <-  rank(mod$clus.avg.silwidths,  ties.method = "first")# large better 
WD <-  rank(-1*mod$within.dist, ties.method = "first") ### mean within distance  - small better  
gap <- rank(-1*(mod$cwidegap), ties.method = "first")  ### with cluster gap .. small better 
SSW <- rank(-1*(mod$within.cluster.ss), ties.method = "first") ### within cluster variance ... small beter 
tb <- cbind(DM, SW,WD,gap,SSW) 
apply(tb, 1, getMode)
}




