#' @title cluster connectivity measure  
#
#' @description Calculates the connectivity validation measure for data points or dissimilarity matrix and
#'  cluster partitioning. Currently uses R's \code{link{dist}} method.  
#  
#
#' @export
clusterConnectivity <- function(distance=NULL, clusters, Data=NULL, neighbSize=10, method="euclidean"){

  if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  nearest <- apply(distance,2,function(x) sort(x,ind=TRUE)$ix[2:(neighbSize+1)])
  nr <- nrow(nearest);nc <- ncol(nearest)
  same <- matrix(clusters,nrow=nr,ncol=nc,byrow=TRUE)!=matrix(clusters[nearest],nrow=nr,ncol=nc)
  conn <- sum(same*matrix(1/1:neighbSize,nrow=nr,ncol=nc))
  return(conn)
}

