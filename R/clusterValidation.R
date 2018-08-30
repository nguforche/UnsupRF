#' @title Cluster Validation by Successive Omission  
#
#' @description Validate a clustering scheme by succerssively omitting the most compact/connected cluster member and  
#'  re-cluster the remaining data points, and continue the process all over.   
#' The hypothesis is that, if the data can be clustered into k optimal cluster, then after
#'  dropping the most compact/connected cluster member, the remaining data points can best be clustered 
#' into k-1 clusters and so on, untill k = 2. By trying different k's, the successive omission algorithm may
#' identify the best k clusters for the data.  
#' 
#' @details
#' This function is experimental and not meant to be used directly by the user
#' Work in progress to validate the procedure 

#
#' @export
successiveOmitClusterValidation  <- function(data, k =  2, vars, cluster.method = "UnsupRF", 
                                             dist.method="euclidean", control=NULL, ensemble = TRUE, sep.measure = "SW" ){
 if(is.null(control)) stop("please provide control paramerter list")
                                             
    Res <- Cluster(data=data[, vars], k=k, cluster.method=cluster.method, control = control)
    clusters <- Res$clusters 
    distance <- Res$distance 
    
    ### get most compact cluster and drop from data 
    ix <- clusterCompactness(distance=distance, clusters = clusters, method = dist.method, ensemble = control$ensemble, 
    compact.measure = control$compact.measure)
    id <-  data[clusters == ix, "id"]
    data <- data[clusters != ix, ]

    ### if k > 2, then re-cluster reduced data into k -1 clusters 
    if(k > 2) {
      sep.mat <- lapply(2:(k-1), function(l) {  
        res <- Cluster(data=data[, vars], k=l, cluster.method=cluster.method, control = control)
        cls <- res$clusters 
        dis <- res$distance 
        ### apply cluster separation to get cluster speration measures 
        clusterSeparation(distance=dis, Data = data[, vars], clusters=cls, method="euclidean")
      })
      ## combine all separation measures and slect the optimal number of clusters 
      ## select by ensemble 
      
      tab  <- do.call(rbind, sep.mat)
      rownames(tab) <- 2:(k-1) 
      ll <- apply(tab, 2, which.max)

      if(ensemble) 
      lopt <- getMode(ll) + 1
      else 
      lopt <- ll[sep.measure] + 1 
      
      if(lopt <= 2)
      return(append(list(id), successiveOmitClusterValidation(data=data, k=2, vars=vars, cluster.method=cluster.method, 
              dist.method=dist.method, control = control, ensemble = ensemble, sep.measure = sep.measure )))
      else {
      return(append(list(id), successiveOmitClusterValidation(data=data, k=lopt, vars=vars, cluster.method=cluster.method, 
               dist.method=dist.method, control = control, ensemble = ensemble, sep.measure = sep.measure )))
      }
    } else {
    ## cluster data 
    res <- Cluster(data=data[, vars], k=2, cluster.method=cluster.method, control = control)
    cls <- res$clusters 
    dis <- res$distance 
    mod <- clusterStats(d = as.matrix(dis), clustering = cls, alt.clustering=NULL)
    sep <- max(mod$max.separation)
    dia <- mod$diameter
    ### if sep > dia then  there are two distinct clusters else single cluster 
    ### how do we determine indicies of the two clusters: just split id by clusters and return one split 
    if(sep > sum(dia)) id <- split(data$id, f = cls)[[1]]
    else id = data$id  
    
    return(list(id))
    }
    
 }
  
  
#' @export
successiveOmitClusterValidationWithDropOut  <- function(data, k =  2, drop = 1, vars, cluster.method = "UnsupRF", 
                                             dist.method="euclidean", control=NULL, ensemble = TRUE, sep.measure = "SW" ){
 if(is.null(control)) stop("please provide control paramerter list")
                                             
    Res <- Cluster(data=data[, vars], k=k, cluster.method=cluster.method, control = control)
    clusters <- Res$clusters 
    distance <- Res$distance 
    
    ### get most compact cluster and drop from data     
    index <- clusterCompactness.ind(distance=distance, clusters = clusters, method = dist.method)

#print(index)

    ### if k > 2, then re-cluster reduced data into k -1 clusters 
    if(k > 2) {    
    ix <- index[1:drop]  
## more than one cluster to drop  and  number of clusters is at least 2 times larger than number of clusters to  drop
## else (not sure the exact outcome here )...  TODO    perhaps break recursion with  return(id) ? what about the case k-drop == 1 ?  
    if(length(ix) > 1){  ### this implies drop > 1
    id <-  lapply(ix, function(xx) data[clusters == xx, "id"]) 
    names(id) <- NULL  
    data <- data[!clusters%in%ix, ]
    if(abs(k-drop)==1) return(id)
    } else {
    id <- data[clusters == ix, "id"]
    data <- data[clusters != ix, ]
   }
   
      sep.mat <- lapply(2:(k-1), function(l) {  
        res <- Cluster(data=data[, vars], k=l, cluster.method=cluster.method, control = control)
        cls <- res$clusters 
        dis <- res$distance 
        ### apply cluster separation to get cluster sparation measures 
        clusterSeparation(distance=dis, Data = data[, vars], clusters=cls, method="euclidean")
      })
      ## combine all separation measures and select the optimal number of clusters 
      tab  <- do.call(rbind, sep.mat)
      rownames(tab) <- 2:(k-1) 
      ll <- apply(tab, 2, which.max)
      
      if(ensemble) 
      lopt <- getMode(ll) + 1
      else 
      lopt <- ll[sep.measure] + 1 
      
      if(lopt <= 2)
      return(append(id, successiveOmitClusterValidationWithDropOut(data=data, k=2, drop = drop, vars=vars, cluster.method=cluster.method, 
              dist.method=dist.method, control = control, ensemble = ensemble, sep.measure = sep.measure)))
      else {
      return(append(id, successiveOmitClusterValidationWithDropOut(data=data, k=lopt, drop = drop,vars=vars, cluster.method=cluster.method, 
               dist.method=dist.method, control = control, ensemble = ensemble, sep.measure = sep.measure)))
      }
      
    } else {
    
    ix <- index[1]      
    id <- data[clusters == ix, "id"]
    data <- data[clusters != ix, ]
   
 ## cluster data 
    res <- Cluster(data=data[, vars], k=2, cluster.method=cluster.method, control = control)
    cls <- res$clusters      
    dis <- res$distance 
    mod <- clusterStats(d = as.matrix(dis), clustering = cls, alt.clustering=NULL)
    sep <- max(mod$max.separation)
    dia <- mod$diameter 
    ### if sep > dia then  there are two distinct clusters else single cluster 
    ### how do we determine indicies of the two clusters: just split id by clusters and return one split 
    if(sep > sum(dia)) id <- split(data$id, f = cls)[[1]] 
    else id = data$id 
   
    return(list(id))
    }
    
 }
  
  
  
  
  
  
  
  
  
  
  
