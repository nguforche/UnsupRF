#' @title Best number of cluster partitions.  
#
#' @description Get "best" number of clustering partitions determined by \code{\link{successiveOmitClusterValidation}} for 
#' k = 2, ..., K using cluster connectivity. Specifically, this function runs \code{\link{successiveOmitClusterValidation}}  
#' for each cluster number in krange and selects the optimal number of cluster by looking at the 
#' clustering with the best separation or compactnest  
#'
#' @param data data matrix to cluster 
#' @param krange integer vector. Numbers of clusters to try
#' @param  dist.method  Use \code{\link{RFdist}} or R's base \code{\link{dist}} method to compute dissimilarity 
#' matrix. Default is "UnsupRF". 
#' @param parallel  run in parallel ? Currently uses \code{\link{mcapply}} which works only on linux operating systems
#' @param mc.cores  number of CPU cores 
#' @param ensemble  take the ensemble of the seperation measures ? 
#' @param control  list with control parameters: 
#' \itemize{
#' \item ntree number of trees for \code{\link{RFdist}}
#' \item no.rep number of repetitions or forests for \code{\link{RFdist}}
#' \item neighbSize  number of nearest neighbors to compute connectivity validation measure. See the function 
#'  \code{\link{clusterConnectivity}}
#' \item method \code{hclust} method 
#' \item RF.parallel  type of parallization for \code{\link{RFdist}}
#' \item compact.measure  one of five compactness measure. See \code{\link{clusterCompactness}}
#' \item combined: if ensemble = FALSE, should the seperation measures be combined?  
#' }
#' @return clustering for data 
#' @details
#' This function is experimental .. work in progress to validate the procedure 
#' @export
bestPartitions  <- function(data, krange=2:10,  cluster.method = "UnsupRF", dist.method="euclidean", 
                        parallel = FALSE, mc.cores = 2, control = NULL, sep.measure = "SW", ensemble = TRUE){
if(is.null(control)) stop("please provide control paramerter list")
                       
data <- as.data.frame(data) 
rhs.vars <- colnames(data)
data[, "id"] <- 1:nrow(data) 

  if(parallel) {
    pfun <-  get("mclapply")
  } else {
    pfun = get("lapply")
  }

### compute dissimilarity 
if(cluster.method == "UnsupRF") {
           distance <- RFdist(data=data[, rhs.vars],  syn.type = "permute", importance= FALSE, parallel = control$RF.parallel, 
           no.rep = control$no.rep, ntree = control$ntree)$RFdist
#                              
} else  { distance <- dist(data[, rhs.vars]) }
id <- data[, "id"]

Res <- pfun(krange, function(xx, ...){

  ss <- successiveOmitClusterValidation(data = data, k = xx, vars = rhs.vars, cluster.method = cluster.method, control = control,
   ensemble = ensemble, sep.measure = sep.measure)  

  ids <-ss 
  ids[[length(ids)+1]] <- setdiff(id, unlist(ss))
  clusters <- numeric(nrow(data))
  for(ii in 1:length(ids))
    clusters[ids[[ii]]] <- ii 
    
con <-  clusterSeparation(distance=distance, Data = data[, rhs.vars], clusters=clusters)
com <-  sum(clusterCompactness.values(distance=distance, clusters = clusters, method = dist.method,
            compact.measure = control$compact.measure))  
  
if(!ensemble){
switch(sep.measure, 
connectivity ={ con <- -1*clusterConnectivity(distance = distance, clusters = clusters, neighbSize= control$neighbSize)}, 
SW = {con <- con[1]}, 
dunn = {con <- con[2]}, 
dunn2 = {con <- con[3]}, 
max.sep = {con <- con[4]}, 
wd.ratio = {con <- con[5]},
Ave.Between = {con <- con[6]}, 
ch = {con <- con[7]}, 
sindex = {con <- con[8]},
DB= {con <- con[9]},
G2={con <- con[10]},
G3={con <- con[11]},
entropy={con <- con[12]},
stop("Wrong value for cluster method")
)
} 
list(clusters = clusters, con = con, com = com)
}, mc.cores = mc.cores)

if(ensemble){
 tab  <- do.call(rbind, lapply(Res, function(xx) xx$con))#[, c("SW", "dunn2", "ch", "sindex", "DB")]
 ll <- apply(tab, 2, which.max)
 ix <- getMode(ll)
} else {
con <-  unlist(lapply(Res, function(xx) xx$con))
	if(control$combined){
	com <-  unlist(lapply(Res, function(xx) xx$com))
	ix1 <- rank(con, ties.method = "first")
	ix2 <- rank(com, ties.method = "first")	
             if(length(ix1) == 2) ix <- ix1[1]  ## prefer well seperated two clustering 
             else ix <- ifelse(ix1[1] == ix2[1], ix1[1], ifelse(ix1[2] == ix2[2], ix1[2], 
                        ifelse(ix1[3]==ix2[3], ix1[3], ix1[1]))) 
	} else 
	ix <- which.max(con)
}            
clusters <- lapply(Res, function(xx) xx$clusters)
clusters[[ix]]
}


#' @export
bestPartitionsWithDropOut  <- function(data, krange=2:10,  drop = 1, cluster.method = "UnsupRF", dist.method="euclidean", 
                        parallel = FALSE, mc.cores = 2, control = NULL, sep.measure = "SW", ensemble = TRUE){
if(is.null(control)) stop("please provide control paramerter list")
                       
data <- as.data.frame(data) 
rhs.vars <- colnames(data)
data[, "id"] <- 1:nrow(data) 

  if(parallel) {
    pfun <-  get("mclapply")
  } else {
    pfun = get("lapply")
  }

### compute dissimilarity 
if(cluster.method == "UnsupRF") {
           distance <- RFdist(data=data[, rhs.vars],  syn.type = "permute", importance= FALSE, parallel = control$RF.parallel, 
           no.rep = control$no.rep, ntree = control$ntree)$RFdist
#                              
} else  { distance <- dist(data[, rhs.vars],  method = dist.method) }
id <- data[, "id"]

Res <- pfun(krange, function(xx, ...){

  ss <- successiveOmitClusterValidationWithDropOut(data = data, k = xx, drop = drop, vars = rhs.vars, 
       cluster.method = cluster.method, control = control,ensemble = ensemble, sep.measure = sep.measure)  

  ids <-ss   
  ids[[length(ids)+1]] <- setdiff(id, unlist(ss))
  clusters <- numeric(nrow(data))

  for(ii in 1:length(ids))
    clusters[ids[[ii]]] <- ii 

#print(clusters)
    
con <-  clusterSeparation(distance=distance, Data = data[, rhs.vars], clusters=clusters)

com <-  sum(clusterCompactness.values(distance=distance, clusters = clusters, method = dist.method, 
            compact.measure = control$compact.measure))  

if(!ensemble){
switch(sep.measure,
connectivity ={ con <- -1*clusterConnectivity(distance = distance, clusters = clusters, neighbSize= control$neighbSize)}, 
SW = {con <- con[1]}, 
dunn = {con <- con[2]}, 
dunn2 = {con <- con[3]}, 
max.sep = {con <- con[4]}, 
wd.ratio = {con <- con[5]},
Ave.Between = {con <- con[6]}, 
ch = {con <- con[7]}, 
sindex = {con <- con[8]},
DB= {con <- con[9]},
G2={con <- con[10]},
G3={con <- con[11]},
entropy={con <- con[12]},
stop("Wrong value for cluster method")
)
}
list(clusters = clusters, con = con, com = com)
}, mc.cores = mc.cores)
if(ensemble){
 tab  <- do.call(rbind, lapply(Res, function(xx) xx$con))
 ll <- apply(tab, 2, which.max) 
 ix <- getMode(ll)
} else {
con <- unlist(lapply(Res, function(xx) xx$con))
	if(control$combined){
	com <-  unlist(lapply(Res, function(xx) xx$com))

	
	ix1 <- rank(con, ties.method = "first")
	ix2 <- rank(com, ties.method = "first")	
             if(length(ix1) == 2) ix <- ix1[1]  ## prefer well seperated two clustering 
             else ix <- ifelse(ix1[1] == ix2[1], ix1[1], ifelse(ix1[2] == ix2[2], ix1[2], 
                        ifelse(ix1[3]==ix2[3], ix1[3], ix1[1]))) 
	} else 
	ix <- which.max(con)
}            
clusters <- lapply(Res, function(xx) xx$clusters)
clusters[[ix]]
}

















