#' @title cluster Separation measures  
#
#' @description Calculate 5 between cluster separation measures  # 
#' @export
clusterSeparation <- function(distance=NULL, Data = NULL, clusters, method="euclidean", DB = TRUE, ...) { 

if(DB & is.null(Data) ) stop("please provide Data to compute DB index")
dat <- Data

 if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
 if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
 if (class(distance)=="dist") distance <- as.matrix(distance)

mod <- clusterStats(d = distance, clustering = clusters, alt.clustering=NULL, sepindex = TRUE, ...)

SW <-  mod$avg.silwidth # large better 
dunn <- mod$dunn ## larger better 
dunn2 <- mod$dunn2 ## larger better
max.sep <- max(mod$max.separation) #mod$sindex ## large better 
wb <- -1*mod$wb.ratio   ### average within distance/average between distance ... smaller the  better so we take the negative to make larger the better  
Ave.Between <- mod$average.between  ## average distance between clusters --> large better 
ch <- mod$ch ## large better 
sindex <- mod$sindex 
DB <- ifelse(DB,  -1*index.DB(x = dat, cl = clusters, d = distance)$DB, NA) ## Davies-Bouldin  smaller better ---> so we take nagative to make largger the better  
G2 <- index.G2(d=distance, cl = clusters) ## adaptation of Goodman and Kruskal's gamma statistics large better
G3 <- -1*index.G3(d=distance, cl = clusters) ## c-index or g3 Hubert and Levin  smaller better ---> note we take the nagative to make larger the better 
entropy <- -1*mod$entropy ### entropy ... smaller the better so we take the nagative  
cbind(SW=SW,dunn=dunn, dunn2=dunn2, max.sep=max.sep, wb.ratio=wb, Ave.Between=Ave.Between, 
ch=ch, sindex=sindex, DB=DB, G2=G2, G3=G3, entropy=entropy) 
}


