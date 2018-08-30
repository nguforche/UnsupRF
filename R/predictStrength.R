#' @title  cluster validation by prediction strength  
#
#' @description Calculate the mediods of a clustering by finding the point that has the minimum average/sum 
#' distance to all other points in the cluster. Proximity between data points can be 
#' provided by \code{RFdist} or using the Gower's general similarity coefficient. 
#
#' @name clvpredictStrength
#
#' @param dat a data matrix    
#' @param dist.mat dissimilarity matrix obtained from the data matrix "dat"    
#' @param method  \code{\link{hclust}} clustering method
#' @param nBoots  number of bootstraps 
#' @param krange integer vector. Numbers of clusters which are to be tried 
#' @param balance  perform balance bootstrap ? 
#' @param parallel run in parallel ?
#' @param mc.cores number of CPU cores 
#' @param OOB  use out-of-bag from bootstrap as test data?  
#' @param seed random seed 
#' @param \dots further arguments passed to or from other methods.
#' @details
#' This function is experimental .. more details coming   
NULL
#' @rdname clvpredictStrength
#' @export
clvpredictStrength <- function(dat, ...) UseMethod("clvpredictStrength")
#' @rdname clvpredictStrength
#' @export
clvpredictStrength.default <- function(dat, dist.mat, method = "ward.D2", nBoots = 20, krange = 2:5, 
                             balanced = FALSE, parallel = FALSE, mc.cores=2, OOB = TRUE, seed = 12345, ...){

if(!inherits(dist.mat, "dist")) stop("Please provide a dissimilarity matrix") 

  if(parallel) {
    set.seed(seed, "L'Ecuyer") 
    pfun <-  get("mclapply")
  } else {
    set.seed(seed)
    pfun = get("lapply")
  }
 
perf.oob <- NULL 
oob <- NULL 
HC.full <- hclust(d= dist.mat, method= method, members= NULL)


n <- attr(dist.mat, "Size")
ix.boot <- 1:n
BOOT <- createResample(ix.boot, nBoots, FALSE)  ## create bootstrap samples 
if(balanced) Boot <- sample(matrix(t(replicate(nBoots, ix.boot)))) 

### mclapply is not available on windows   
Boot.res <-  pfun(1:nBoots, function(kk, ...){

inbag <- BOOT[, kk]
if(balanced) inbag <- BOOT[((kk-1)*n+1):(kk*n)]
dat.trn <- dat[inbag, ,drop = FALSE]
distmat <- as.matrix(dist.mat)
d.trn <- as.dist(distmat[inbag, inbag])
hclust.trn <- hclust(d=d.trn, method= method, members= NULL)

if(OOB){
outbag <- setdiff(ix.boot, inbag)
dat.tst <- dat[outbag, ,drop = FALSE]
d.tst <- as.dist(distmat[outbag, outbag])
hclust.tst <- hclust(d=d.tst, method= method, members= NULL)
}

##### 
clus.perf <- lapply(krange, function(xx){
full.clusters <- factor(cutree(HC.full, xx))
trn.clusters <- cutree(hclust.trn, xx)

centers <- mediod(x = d.trn, clusters=trn.clusters)
centers$mediod <- dat.trn[centers$ids, ,drop = FALSE]		

l1 <- factor(predict(object=centers, org.data=NULL, newdata = dat))
#### performance 
dd <- data.frame(obs = full.clusters, pred = l1, stringsAsFactors = TRUE) 
levels(dd$pred) <- levels(full.clusters)
perf <- multiClassSummary(data = dd, lev = levels(full.clusters))

if(OOB){
tst.clusters <- factor(cutree(hclust.tst, xx))
l2 <- factor(predict(object=centers, org.data=NULL, newdata = dat.tst))
d.oob <- data.frame(obs = tst.clusters, pred = l2, stringsAsFactors = TRUE) 
levels(d.oob$pred) <- levels(tst.clusters)
perf.oob <- multiClassSummary(data = d.oob, lev = levels(tst.clusters))
perf.oob <- cbind(k = xx, t(perf.oob))
}

list(full = cbind(k = xx, t(perf)), oob = perf.oob)
})
full <- do.call(rbind, lapply(clus.perf, function(ss) ss$full))
if(OOB) oob <-do.call(rbind, lapply(clus.perf, function(ss) ss$oob)) 
list(full = full, oob = oob)
})
alpha = 0.05 
Tab <- data.frame(do.call(rbind, lapply(Boot.res, function(ss) ss$full)))
Mn <- ddply(Tab, .variables = "k", .fun = numcolwise(mean) , na.rm = TRUE)
Sd <- ddply(Tab, .variables = "k", .fun = numcolwise(mean) , na.rm = TRUE)
CI <- ddply(Tab, .variables = "k", numcolwise(quantile), probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
Full <- list(Mean = Mn, Sd = Sd, CI = CI) 
if(OOB){
Tab <- data.frame(do.call(rbind, lapply(Boot.res, function(ss) ss$oob)))
Mn <- ddply(Tab, .variables = "k", .fun = numcolwise(mean) , na.rm = TRUE)
Sd <- ddply(Tab, .variables = "k", .fun = numcolwise(mean) , na.rm = TRUE)
CI <- ddply(Tab, .variables = "k", numcolwise(quantile), probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
oob <- list(Mean = Mn, Sd = Sd, CI = CI) 
}
res <- list(Full = Full, OOB = oob)
class(res) <- "clvpredictStrength"
res
}


#' @rdname  clvpredictStrength
#' @method plot clvpredictStrength
#' @export

plot.clvpredictStrength  <- function(x, perf.measure = "Sensitivity") {
  if (!inherits(x, "clvpredictStrength")) stop("Object must be a \"clvpredictStrength \"'")
means <- x$Full$Mean
Sd <- x$Full$Sd
tab <- data.frame(k = means$k, value = means[, perf.measure], ymin = means[, perf.measure] - Sd[,perf.measure], 
                 ymax = means[, perf.measure] + Sd[,perf.measure])

pp <- ggplot( ) +
  geom_line(data = tab, aes(x = k , y = value), colour = "blue", size = 1.3) + 
  geom_point(data = tab, aes(x = k, y = value),  colour = "blue", position = position_dodge(width = 0.02), size = 3) + 
  geom_errorbar(data = tab, aes(x = k, ymin = ymin, ymax = ymax), position = position_dodge(width = 0.3), width = 0.2) + 
#  geom_vline(xintercept = 5, colour = "darkgreen") + 
  scale_x_continuous(name="Number of clusters",breaks=2:30) + ylab(perf.measure) +  
  labs(title = "Estimating the number of clusters via the prediction strength") + 
  theme(axis.title.x=element_text(size=14,face="bold"), 
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(size=14,face="bold"), 
        axis.text.x = element_text(size = 13, face="bold",colour = "gray40"),
        legend.title = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 13, face="bold",colour = "gray40")) 
print(pp)
}


#' @rdname  clvpredictStrength
#' @method print clvpredictStrength
#' @export
print.clvpredictStrength <- function(x,...){
  if (!inherits(x, "mediod")) stop("Object must be a \"mediod \"'")
  print("*** Cluster mediod row indices in original data ***")
  print(x$ids)
}


