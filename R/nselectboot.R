 #' Selection of optimal number of clusters using bootstrap.
 #
 #' This is a minor modification of the  nselectboot function in the package fpc, so you can use package::nselectboot
 #' to run the perferred version if both packaes are loaded. See \code{\link[fpc]{nselectboot}} for details.
 #
 #' @param data data or dissimilarity matrix
 #' @param B number of boostrap
 #' @param distances (logical) is data a dissimilarity matrix ?
 #' @param clustermethod clustering method, see \code{\link[fpc]{nselectboot}}
 #' @param classification see \code{\link[fpc]{nselectboot}}
 #' @param  krange integer vector; numbers of clusters to be tried
 #' @param count logical. If TRUE, numbers of clusters and bootstrap runs are printed.
 #' @param nnk number of closest neighbors for classification
 #' @param \dots further arguments passed to or from other methods.
#  
#' @export
nselectboot <- function (data, B = 50, distances = inherits(data, "dist"), clustermethod = NULL, 
                           classification = "averagedist", krange = 2:10, count = FALSE, 
                           nnk = 1, ...) 
{
  dista <- distances
  data <- as.matrix(data)
  if (classification == "average") {
    if (dista) 
      dmat <- data
    else dmat <- as.matrix(dist(data))
  }
  stab <- matrix(0, nrow = B, ncol = max(krange))
  n <- nrow(data)
  for (k in krange) {
    if (count) 
      cat(k, " clusters\n")
    for (i in 1:B) {
      if (count) 
        print(i)
      d1 <- sample(n, n, replace = TRUE)
      d2 <- sample(n, n, replace = TRUE)
      if (dista) {
        dmat1 <- data[d1, d1]
        dmat2 <- data[d2, d2]
      }
      else {
        dmat1 <- data[d1, ]
        dmat2 <- data[d2, ]
      }
      clm1 <- clustermethod(dmat1, k = k, ...)
      clm2 <- clustermethod(dmat2, k = k, ...)
      cj1 <- cj2 <- rep(-1, n)
      cj1[d1] <- clm1$partition
      cj2[d2] <- clm2$partition
      if (dista) {
        
        if (identical(clustermethod, pamkCBI)) {
          cj1 <- classifdist(data, cj1, method = classification, 
                             centroids = clm1$result$pamobject$medoids, nnk = nnk)
          cj2 <- classifdist(data, cj2, method = classification, 
                             centroids = clm2$result$pamobject$medoids, nnk = nnk)
        }
        
        if (identical(clustermethod, claraCBI)) {
          cj1 <- classifdist(data, cj1, method = classification, 
                             centroids = clm1$result$medoids, nnk = nnk)
          cj2 <- classifdist(data, cj2, method = classification, 
                             centroids = clm2$result$medoids, nnk = nnk)
        }
        
        if (identical(clustermethod, kmeansCBI)) {
          centroids1 <- clm1$result$centers
          centroids2 <- clm2$result$centers
          cj1 <- classifnp(data, cj1, method = classification, 
                           centroids = centroids1, nnk = nnk)
          cj2 <- classifnp(data, cj2, method = classification, 
                           centroids = centroids2, nnk = nnk)                   
        }
        
      }
      else {
        centroids <- NULL
        if (classification == "centroid") {
          if (identical(clustermethod, kmeansCBI)) {
            centroids1 <- clm1$result$centers
            centroids2 <- clm2$result$centers
          }
          if (identical(clustermethod, claraCBI)) {
            centroids1 <- clm1$result$pamobject$medoids
            centroids2 <- clm2$result$pamobject$medoids
          }
          if (identical(clustermethod, pamkCBI)) {
            centroids1 <- clm1$result$pamobject$medoids
            centroids2 <- clm2$result$pamobject$medoids
          }
        }
        cj1 <- classifnp(data, cj1, method = classification, 
                         centroids = centroids1, nnk = nnk)
        cj2 <- classifnp(data, cj2, method = classification, 
                         centroids = centroids2, nnk = nnk)
      }
      ctable <- table(cj1, cj2)
      nck1 <- rowSums(ctable)
      stab[i, k] <- sum(nck1^2 - rowSums(ctable^2))
    }
  }
  stab <- stab/n^2
  stabk <- rep(NA, max(krange))
  for (k in krange) stabk[k] <- mean(stab[, k])
  kopt <- which.min(stabk)
  out <- list(kopt = kopt, stabk = stabk, stab = stab)
}

