#' The Jaccard similarity coefficient is defined as:
#' \deqn{J = \frac{n_{11}}{n_{11} + n_{10} + n_{01}}}.
#'
#' In the special case that the Jaccard coefficient results in \eqn{0/0},
#' we define \eqn{J = 0}. For instance, this case can occur when both clusterings
#' consist of all singleton clusters.
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.

#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Jaccard coefficient for the two sets of cluster labels (See' Details.)
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Jaccard similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' jaccardIndex(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Jaccard similarity coefficient between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' jaccardIndex(iris_kmeans, iris_hclust)
#' }
jaccardIndex <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  jaccard_out <- with(com_table, n_11 / (n_11 + n_10 + n_01))

  # In the case where 'labels1' and 'labels2' contain all singletons, the Jaccard
  # coefficient results in the expression 0 / 0, which yields a NaN value in R.
  # We define such cases as 0.
  if (is.nan(jaccard_out)) {
    warning("The two clusterings contain all singletons -- returning 0.")
    jaccard_out <- 0
  }
  jaccard_out
}


