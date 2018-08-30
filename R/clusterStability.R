#' Cluster Stability. 
#' \deqn{J = \frac{n_{11}}{n_{11} + n_{10} + n_{01}}}.
#'
#' Computes the cluster stability measure defined in Famili et.al based on the 
#' idea of cluster immovability on partition. Cluster immovability is the rate at which 
#' the content of a cluster remains unchanged during the clustering process. 
#'
#' @name clusterStability

#' @param c1 a vector of \code{n} clustering labels
#' @param c2 a vector of \code{n} clustering labels
#' @return the Jaccard coefficient for the two sets of cluster labels (See
#' Details.)
#' @references
#' Famili, A. Fazel, Ganming Liu, and Ziying Liu. "Evaluation and optimization of clustering in 
#' gene expression data analysis." Bioinformatics 20.10 (2004): 1535-1545.
#
NULL
#' @rdname clusterStability
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Jaccard similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' jaccard_indep(labels1, labels2)
#' 
#' # Here, we cluster the iris data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Jaccard similarity coefficient between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' jaccard_indep(iris_kmeans, iris_hclust)
#' }
#
#' @export
Stability <- function(c1, c2) {
length(intersect(c1,c2))/length(c1)
}


clusterStability <- function(objects, Hmod, krange=2:5, k = NULL ) {
n <- max(krange) 
## c in krange 
# k  in 1, 2, ..., n-c 
# i = (c+1):(c+k) 

  
## now do l partitions
GS.c <- sapply(krange, function(cc){ 
P.c = cutree(Hmod, k = cc)
Mems.c <- split(objects, f = P.c) 

if(is.null(k)) kk = n-cc
else {
if(k <= n-cc) kk = k
else kk = n-cc 
} 

min.j <- sapply(unique(P.c), function(l){
max.i <- sapply((cc+1):(cc+kk), function(i){   
P.i = cutree(Hmod, k = i)
Mems.i <- split(objects, f = P.i) 
max(sapply( unique(P.i), function(j) Stability (Mems.c[[l]], Mems.i[[j]]) ))
}
)
max.i <- max.i/length(Mems.c[[l]])
min(max.i)  ## min of k thresholds for each partition c and membership l 
}
)
sum(min.j)/cc
}
)
GS.c
}



















