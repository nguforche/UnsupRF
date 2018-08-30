#' @title  mediods of a clustering procedure
#
#' @description Calculate the mediods of a clustering by finding the point that has the minimum average/sum 
#' distance to all other points in the cluster. Proximity between data points can be 
#' provided by \code{RFdist} or using the Gower's general similarity coefficient. 
#
#' @name mediod 
#
#' @param x a dissimilarity matrix, or a data frame/ matrix   
#' @param clusters clustering 
#' @param data data frame/ matrix   
#' @param fun character name of the function to determine minimum 
#' distance between points in a cluster. Should be either \code{mean}, 
#' \code{median} or \code{sum} (default)
#' @param weights ptional vector of weights for variables in data. See \code{growdis} in the 
#' \code{FD} package 
#' @param \dots further arguments passed to or from other methods.
#' @details
#' \code{mediod} is the main function to compute the mediod, while 
#' \code{mix.dist} computes the distance between obersavations based on veriables 
#' of mix-types: binary, categorical, and continuous using the Gower's general similarity 
#' coefficient.  
#' @return data matrix of cluster mediods and the corresponding row indicies of 
#' the mediods in the original data
#' @references
#' Gower, John C. "A general coefficient of similarity and some of its properties." Biometrics (1971): 857-871.
#
#' @importFrom FD gowdis
NULL
#' @rdname mediod
#' @export
mediod <- function(x, ...) UseMethod("mediod")
#' @rdname mediod
#' @export
mediod.default <- function(x, clusters, fun = "sum", weights, ...){
centers <- NULL  

if(inherits(x, "dist")){ 
data <- as.matrix(x)
} else {
 data <- as.matrix(mix.dist(x, weights))
}
# split dissimilarity matrix by clusters: to produce sub-dissimilarity matrices for each 
# cluster and  find data point with minimum average/median/sum of dissimilarity to 
# every other data point in cluster 
  
  clusters <- factor(clusters) 
  levels(clusters) <- paste0("cluster:", levels(clusters))
  nme  <- levels(clusters)
  ix <- split(1:nrow(data), f = clusters)
  X <- split.data.frame(data, f = clusters)
  
centers.ids <- sapply(nme, function(x){
    tb <- X[[x]][, ix[[x]]]
    m <- colnames(tb)
    iy <- which.min(apply(tb, 2, get(fun) ))
    m[iy]
  }) 

if(!inherits(x, "dist")) centers <- x[centers.ids, ]
 
res <- list(mediod = centers, ids = centers.ids)    
class(res) <- "mediod"
return(res)
}
#
#' @rdname  mediod 
#' @export
mix.dist <- function(data, weights){
  dc <- sapply(data, data.class)
  if (all(dc == "numeric")) 
    D <- dist(data)
  else {
    if (any(dc == "logical")) 
      data[, dc == "logical"] <- sapply(data[, dc == "logical"], as.numeric)
    fac <- sapply(data[, dc == "factor", drop = FALSE], function(x) length(levels(x)))
    bin <- names(fac)[fac == 2]
    data[, bin] <- sapply(data[, bin], function(x) as.numeric(x) - 1)
    D <- gowdis(x = data, w = weights, ord = "metric")
    D <- sqrt(D)
} 
return(D)
}

#' @rdname  mediod 
#' @method print mediod
#' @export
print.mediod <- function(x,...){
  if (!inherits(x, "mediod")) stop("Object must be a \"mediod \"'")
  print("*** Cluster mediod row indices in original data ***")
  print(x$ids)
}



  
