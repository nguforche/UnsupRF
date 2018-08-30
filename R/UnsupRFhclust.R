#' @title Unsupervised random forest with hclust: cluster data and evaluate predictive strength of clusters.
#
#' @description This function takes a dissimilarity matrix, such as the Random Forest dissimilarity matrix 
#' from \code{\link{RFdist}} and contructs a hirearchical clustering object using the \code{\link{hlust}} package. 
#'   It then evaluates the predictive ability of different clusterings k = 2:K by 
#'  predicting a binary response variable based on cluster memberships. The results can be used to 
#' validate and select the best number of clusters.  See Ngufor et al. 
#' 
#'  It takes a standard  \code{formula}, a data matrix \code{dat} containing the binary response, and a 
#' disimilarity  matrix \code{rfdist} derived from \code{dat}  and computes the AUCs for 
#' three logistic regression models: 
#'     (1) a model with the provided predictors in the formula, (2) a model with 
#' k = 1:K clusters generated from the \code{hclust} object, 
#'     and (3)  a model with k=1:K randomly generated clusters. 
#'     This is then repeated over \code{cv} cross-validation. See example below. 
#' 
#' Note that \code{dat} must be used to compute the dissimilarity matrix and passed to the function unchanged. 
#'   Otherwise there is no guarantee that the method will accurately assign clusters to the right observations 
#' in \code{dat}. Perhaps there might be a way to check for this?  
#'   
#' @name UnsupRFhclust  
# 
#' @param formula an R formuala. Note, only binary outcomes are currently supported. 
#' @param dat data matrix 
#' @param rfdist  dissimilarity matrix, such as the RF distance matrix computed using \code{\link{RFdist}} 
#'     based on the \code{dat} data matrix.
##' @param Hclust.model hierarchical clustering model  
#' @param hclust.method the agglomeration method to be use see \code{link{hlust}}. 
#' @param  RFdist  RF distance matrix computed from \code{\link{RFdist}}.
#' @param K maximum number of clusters to generate for hierarchical clustering model 
#' @param cv number of cross-validation to perform using data, must be at least 2 
#' @param parallel (logical) run in parallel ? 
#' @param mc.cores number of cores 
#' @param seed random seed  
#' @param caret.method  classifcation  method to use in caret: "glm" and "rf" currently tested. 
#' @param fitControl  Control the computational nuances of the caret \code{\link{train}} function. 
#'  See the \code{\link{caret}} package 
#' @param \dots further arguments passed to or from other methods.
#' @return  a data frame with columns: 
#' \enumerate{
#' \item UnsupRFhclust and UnsupRFhclust.caret
#' \enumerate{
#'   \item Hclust: hierarchical clustering model 
#'   \item perf: a data frame with columns: AUC, cluster (cluster numbers), CV (cross-validation number), and type 
#' (for one of the three types of models mention in the description)  
#'  }
#' \item cluster_internal_validation
#' \enumerate{
#'   \item clusters: cluster memberships for the optimal number of clusters (out put of cluster_internal_validation)
#'   \item kopt: best number of clusters obtained by majprity rule on the table results below. 
#'   \item table: matrix of internal validation metrics  (out put of cluster_internal_validation) with columns 
#'         \itemize{
#'         \item sep: cluster seperation. Higher values indicates better clustering 
#'         \item toother: sum of the average distances of a point in a cluster to points in other clusters. Higher values 
#'               indicates better clustering  
#'         \item within : sum of average distance within clusters: Smaller values indicates better clustering  
#'         \item between: average distance between clusters. Higher values indicates better clustering 
#'         \item ss: within clusters sum of squares. Smaller values indicates better clustering 
#'         \item silwidth: average silhouette width. Higher values indicates better clustering. See \code{\link{silhouette}}. 
#'         \item dunn: Dunn index. Higher values indicates better clustering
#'         \item dunn2: Another version of Dunn index
#'         \item wb.ratio: (negative) ratio of  average distance between clusters to average distance within clusters.  
#'                Smaller values indicates better clustering (or positive and large value is better)  
#'        \item ch: Calinski and Harabasz index. Higher values is better 
#'        \item entropy: negative entropy. Smaller values are better 
#'        \item w.gap:  sum of vector of widest within-cluster gaps. Small is better 
#'         }
#'   }
#' }
#' @details
#' \enumerate{
#'  \item \code{UnsupRFhclust} evaluates the predictive strength of the clusters using base \code{glm} 
#'  directly while \code{UnsupRFhclust.caret} uses glm through the caret package (caret package required). This offers 
#' an interphase to use different classifcation models in the caret package.  
#'  \item  The function \code{cluster_internal_validation} takes a hierarchical clustering model and a 
#'    dissimilarity matrix (e.g output of \code{\link{RFdist}}, and runs the cluster.stats function in the 
#' \code{fpc} package for K  different number of clusters and compute several 
#' internal cluster validation metrics and selects the best number of clusters by majority rule 
#' } 
#' 
#' @references
#' Ngufor, C., Warner, M. A., Murphree, D. H., Liu, H., Carter, R., Storlie, C. B., & Kor, D. J. (2017). 
#' Identification of Clinically Meaningful Plasma Transfusion Subgroups Using Unsupervised Random Forest Clustering. 
#' In AMIA Annual Symposium Proceedings (Vol. 2017, p. 1332). American Medical Informatics Association.
#' 
NULL 
#' @rdname UnsupRFhclust  
#' @export
#' @examples
#' \dontrun{
#' require(plyr)
#' require(ggplot2)
#' data(iris)
#' dat <- iris
#'  # get Random forest dissimilarity matrix 
#' RF.dist <- RFdist(data=dat[, -5], ntree = 10, no.rep=20, 
#'            syn.type = "permute", importance= FALSE)
#
#' form <- as.formula(paste0("Species ~ ", 
#'                paste0(setdiff(names(dat),c("Species")),collapse = "+")))
#' 
#'  # UnsupRFhclust 
#' res <- UnsupRFhclust(formula=form, dat=dat, rfdist = RF.dist$RFdist, K =20,  
#'                      parallel = FALSE, cv = 5, seed = 123)
#' tb <- ddply(res$perf, .variables = c("cluster", "type"), .fun = numcolwise(mean) )
# 
#'  pp <- ggplot( ) +
#'  geom_line(data = tb, 
#'  aes(x = cluster , y = AUC, colour = type, linetype = type), size = 1.3) + 
#'  scale_color_manual(values= c("darkgreen", "darkred", "blue")) + 
#'  geom_vline(xintercept = 3, colour = "darkgreen") + 
#'  scale_x_continuous(name="Number of clusters",breaks=2:30) + ylab("AUC") +  
#'  ylim(c(0.55, 1)) +  
#'  theme(axis.title.x=element_text(size=14,face="bold"), 
#'  axis.title.y=element_text(size=14,face="bold"),
#'  legend.text = element_text(size=14,face="bold"), 
#'  axis.text.x = element_text(size = 13, face="bold",colour = "gray40"),
#'  legend.title = element_text(size=14,face="bold"),
#'  axis.text.y = element_text(size = 13, face="bold",colour = "gray40")) + 
#'  scale_linetype_manual(values=c("solid", "solid", "dotted"))
#' print(pp)
#' 
#' 3 appears to be the best number of clusters 
#' clusters <- cutree(res$Hclust, 3)  
#' 
#' # cluster_internal_validation example  
#' # get hclust object from  UnsupRFhclust 
#'   HCmod = res$Hclust 
#'  rr <- cluster_internal_validation(K = 20, Hclust.model = HCmod,  
#'   RFdist = RF.dist$RFdist, seed = 1234)
#'  rr$table  
#'  rr$kopt 
#' }


### train binary logistic regression model for different clustering from a hierarchical clustering and 
### computes the AUC. 
UnsupRFhclust <- function(formula, dat, rfdist,  hclust.method = "ward.D2", K = 10,  parallel = TRUE, cv = 5, 
                         mc.cores = 2, seed = 12345, verbos = TRUE, ...){

 
  if(parallel) {
    set.seed(seed, "L'Ecuyer") 
    pfun <-  get("mclapply")
  } else {
    set.seed(seed)
    pfun = get("lapply")
  }
  
 ###generate hierarchical cluster object  
Hclust.model  <- hclust(d= rfdist, method=hclust.method)
   
  if(cv <= 1) {cv = 2; warning("cv should be >= 2. Using cv = 2 -fold cross-validations")}   
  xx <- all.vars(formula[[2]])
  ix.cv <- createFolds(dat[, xx], k= cv, list=FALSE)
  dat[, xx] <- factor(dat[, xx])
  if(nlevels(dat[,xx])>2) stop("Only binary classification is currently supported")

  Res <- lapply(1:cv, function(ii){
    
    res <- pfun(2:K, function(kk, ...){
      dat$cluster <- factor(cutree(Hclust.model, kk))  
      ix <-   sample(kk, nrow(dat), replace=TRUE)  
      dat$random.cluster <- factor(ix, levels = unique(ix))
      dat.trn <- dat[ix.cv != ii, ]
      dat.tst <- dat[ix.cv == ii, ] 
     
      mod <- glm(formula=formula, family = binomial, data=dat.trn, control = list(maxit = 100))
      trn <-  predict(mod, newdata = dat.trn, type =  "response")
      tst <-  predict(mod, newdata = dat.tst, type =  "response")
      
      y <- ifelse(dat.trn[, xx] == "Yes", 1, 0) 
      thresh <- as.numeric(PerformanceMeasures(obs = y, pred = as.numeric(trn) )$threshold)
      y <- ifelse(dat.tst[, xx] == "Yes", 1, 0) 
      A1 <- PerformanceMeasures(obs = y, pred = as.numeric(tst), threshold = thresh)$AUC

      form <- as.formula(paste0(paste0(xx, " ~"), paste0("cluster", collapse= "+")))
      mod <- glm(formula=form, family = binomial, data=dat.trn, control = list(maxit = 100))      
      mod$xlevels$cluster <- union(union(mod$xlevels$cluster, levels(dat.trn$random.cluster)), 
              levels(dat.tst$random.cluster))
      trn <-  predict(mod, newdata = dat.trn, type =  "response")
      tst <-  predict(mod, newdata = dat.tst, type =  "response")
            
      y <- ifelse(dat.trn[, xx] == "Yes", 1, 0) 
      thresh <- as.numeric(PerformanceMeasures(obs = y, pred = as.numeric(trn) )$threshold)
      y <- ifelse(dat.tst[, xx] == "Yes", 1, 0) 
      A2 <- PerformanceMeasures(obs = y, pred = as.numeric(tst), threshold = thresh)$AUC
      
      form1 <- as.formula(paste0(paste0(xx, " ~"), paste0("random.cluster", collapse= "+")))
      mod1 <- glm(formula=form1, family = binomial, data=dat.trn, control = list(maxit = 100))
      mod1$xlevels$random.cluster <- union(mod1$xlevels$random.cluster, levels(dat.tst$random.cluster))
      
      tst <-  predict(mod1, newdata = dat.tst, type =  "response")
      thresh <- 0.5
      y <- ifelse(dat.tst[, xx] == "Yes", 1, 0) 
      A3 <- PerformanceMeasures(obs = y, pred = as.numeric(tst), threshold = thresh)$AUC
      
      t1 <- data.frame(CV = ii, cluster = kk, AUC = A1, type = "Predictors")
      t2 <- data.frame(CV = ii, cluster = kk, AUC = A2, type = "H Clusters")
      t3 <- data.frame(CV = ii, cluster = kk, AUC = A3, type = "Random Clusters")
     if(verbos)  cat("done cluster: ", kk, "\n")
      return(rbind(t1, t2, t3))
    }, mc.cores = mc.cores)
    return(do.call(rbind, res))
  })
tab <-   do.call(rbind, Res)
list(Hclust = Hclust.model, perf = tab)
}

#' @rdname UnsupRFhclust  
#' @export
## use caret package 
UnsupRFhclust.caret <- function(formula, dat, rfdist, hclust.method = "ward.D2", K = 10,  parallel = TRUE, cv = 5, 
                        mc.cores = 2, seed = 12345, caret.method = "glm", 
                        fitControl = trainControl(method = "none", classProbs = TRUE), verbos = TRUE, ...){
 print(caret.method)
# suppressMessages(requireNamespace(package="caret", quietly = TRUE))
  
  if(parallel) {
    set.seed(seed, "L'Ecuyer") 
    pfun <-  get("mclapply")
  } else {
    set.seed(seed)
    pfun = get("lapply")
  }
  
  ###generate hierarchical cluster object  
Hclust.model  <- hclust(d= rfdist, method=hclust.method)
 
  if(cv <= 1) {cv = 2; warning("cv should be >= 2. Using cv = 2 -fold cross-validations")}   
  xx <- all.vars(formula[[2]])
  ix.cv <- createFolds(dat[, xx], k= cv, list=FALSE)
  dat[, xx] <- factor(dat[, xx])
  if(nlevels(dat[,xx])>2) stop("Only binary classification is currently supported")

#  fitControl <- caret::trainControl(caret.tranControl)  
     
  Res <- lapply(1:cv, function(ii){
    
    res <- pfun(2:K, function(kk, ...){
      dat$cluster <- factor(cutree(Hclust.model, kk))      
      dat$random.cluster <- factor(sample(kk, nrow(dat), replace=TRUE))
      dat.trn <- dat[ix.cv != ii, ]
      dat.tst <- dat[ix.cv == ii, ] 
 
      mod <- caret::train(formula, data = dat.trn, method = caret.method,  
                   trControl = fitControl, metric = "ROC", ...)                                     
      trn <-  predict(mod, newdata = dat.trn, type =  "prob")[, "Yes"]
      tst <-  predict(mod, newdata = dat.tst, type =  "prob")[,"Yes"]
      
      y <- ifelse(dat.trn[, xx] == "Yes", 1, 0) 
      thresh <- as.numeric(PerformanceMeasures(obs = y, pred = as.numeric(trn) )$threshold)
      y <- ifelse(dat.tst[, xx] == "Yes", 1, 0) 
      A1 <- PerformanceMeasures(obs = y, pred = as.numeric(tst), threshold = thresh)$AUC

      form <- as.formula(paste0(paste0(xx, " ~"), paste0("cluster", collapse= "+")))        
      mod <- caret::train(form, data = dat.trn, method = caret.method,  
                   trControl = fitControl, metric = "ROC", ...)
      trn <-  predict(mod, newdata = dat.trn, type =  "prob")[, "Yes"]
      tst <-  predict(mod, newdata = dat.tst, type =  "prob")[,"Yes"]
      y <- ifelse(dat.trn[, xx] == "Yes", 1, 0) 
      thresh <- as.numeric(PerformanceMeasures(obs = y, pred = as.numeric(trn) )$threshold)
      y <- ifelse(dat.tst[, xx] == "Yes", 1, 0) 
      A2 <- PerformanceMeasures(obs = y, pred = as.numeric(tst), threshold = thresh)$AUC
      
      form <- as.formula(paste0(paste0(xx, " ~"), paste0("random.cluster", collapse= "+")))  
      mod <-  caret::train(form, data = dat.trn, method = caret.method,  
                    trControl = fitControl, metric = "ROC", ...)   
      tst <-  predict(mod, newdata = dat.tst, type =  "prob")[,"Yes"]
      thresh <- 0.5
      y <- ifelse(dat.tst[, xx] == "Yes", 1, 0) 
      A3 <- PerformanceMeasures(obs = y, pred = as.numeric(tst), threshold = thresh)$AUC
      
      t1 <- data.frame(CV = ii, cluster = kk, AUC = A1, type = "Predictors")
      t2 <- data.frame(CV = ii, cluster = kk, AUC = A2, type = "H Clusters")
      t3 <- data.frame(CV = ii, cluster = kk, AUC = A3, type = "Random Clusters")
      if(verbos) cat("done cluster: ", kk, "\n")
      return(rbind(t1, t2, t3))
    }, mc.cores = mc.cores)
    return(do.call(rbind, res))
  })
  tab <-   do.call(rbind, Res)
list(Hclust = Hclust.model, perf = tab)
}
#' @rdname UnsupRFhclust
#' @export
#
#### several cluster internal validation metrics 
cluster_internal_validation <- function(K = 10, Hclust.model,  RFdist, parallel = TRUE,  
                               mc.cores = 2, seed = 1234, ...){

if(parallel) {
    set.seed(seed, "L'Ecuyer") 
     pfun <-  get("mclapply")
  } else {
    set.seed(seed)
     pfun = get("lapply")
  }
  
Res <- pfun(2:K, function(kk, ...){

cluster <- cutree(Hclust.model, kk)
mod <- cluster.stats(d = RFdist, clustering = cluster, alt.clustering=NULL)
sep <-   sum(mod$separation)  ## large better 
toother <- sum(mod$average.toother) ## large better
within <- -1*sum(mod$average.within) ## smaller the better
between <- mod$average.between ## large better
ss <-  -1*mod$within.cluster.ss## smaller the better
silwidth <- mod$avg.silwidth
dunn <- mod$dunn 
dunn2 <- mod$dunn2
#entropy <- mod$entropy 
wb.ratio <- -1*mod$wb.ratio 
ch <- mod$ch ## large better
entropy <- -1*mod$entropy ### large uncertainty not good 
w.gap <- -1*sum(mod$cwidegap)
cbind(sep=sep, toother = toother, within = within, between = between, ss = ss, 
silwidth = silwidth, dunn = dunn, dunn2 = dunn2, wb.ratio = wb.ratio, ch = ch, 
entropy = entropy, w.gap = w.gap)
}, mc.cores = mc.cores)

tab  <- do.call(rbind, Res)
k <- getMode(apply(tab, 2, which.max)+1)
list(clusters = cutree(Hclust.model, k), kopt = k, table = tab) 
}



