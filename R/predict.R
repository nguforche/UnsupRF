#' @title predict cluster memberships for new data. 
#
#' @description Predict cluster memberships for new data by locating points closest to the mediods of a
#' a clustering procedure. Distance between mediod and  newdata is based on the Gower's general similarity 
#' coefficient, which applies to mix-type data.   
# 
#' @param object object of class \code{\link{mediod}}. This requires \code{mediod} to 
#' contain the data matrix of the mediod.
#' @param org.data (optional) the original data used to generate the mediods. 
#' @param newdata A new input data frame.
#' @param ... Further arguments passed to mix.dist.
#' @return clustering for new data  
#
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' data(iris)
#' dat <- iris[, -5]
#' RF.dist <- RFdist(data=dat, ntree = 10, no.rep=20, syn.type = "permute", 
#'                importance=FALSE)
#' # run hclust 
#' HCmod <- hclust(d=RF.dist$RFdist, method="ward.D2", members= NULL)
#' clusters <- cutree(HCmod, 3)   
#' # get mediods of clustering 
#'  med <- mediod(x = RF.dist$RFdist, clusters=clusters, fun = "sum")
#'  print(med)
#' # predict clusters 
#' pred.clusters <- predict(object= med, org.data = dat, newdata = dat) 
#' table(clusters, pred.clusters)
#' }
#
predict.mediod <- function(object, org.data = NULL, newdata,  ...){
if (!inherits(object, "mediod")) stop("Object must be of class  \"mediod \"'")
mediods <- object$mediod
if(is.null(mediods)) {
 if(is.null(org.data)) stop("please provide the original data") 
 else  
 mediods <- org.data[object$ids, ] 
 }
newdata <- data.matrix(newdata)
Mat <- matrix(NA, nrow=nrow(newdata), ncol=nrow(mediods))
newdata <- data.matrix(newdata)
for(ii in 1:nrow(mediods)){
  d <- rbind(newdata, mediods[ii, ])
  d <- as.matrix(mix.dist(d, ...))  
  Mat[, ii] <- d[nrow(d), -ncol(d)] 
}  
apply(Mat, 1, which.min)  
}
















