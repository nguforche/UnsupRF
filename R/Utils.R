optthresh <- function (prob, obs, opt.methods = 9) 
{
  thresh = 0.5
  if (length(unique(obs)) > 1) {
    obs <- as.numeric(as.factor(obs)) - 1
    SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, 
                               Predicted = prob)
    thresh <- optimal.thresholds(SIMDATA, threshold = 101, 
                                 which.model = 1, opt.methods = opt.methods)
    thresh <- ifelse(length(thresh["Predicted"]) >= 1, as.numeric(thresh["Predicted"]), 
                     0.5)
  }
  return(thresh)
}

PerformanceMeasures <- function(obs, pred, threshold=NULL){
  nme = c("AUC", "sensitivity", "specificity")
  if(is.null(threshold)){
    threshold <- optthresh(prob= pred, obs = obs)
  }
  xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)                     
  AUC <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, "AUC"]
  return(cbind.data.frame(AUC = AUC, threshold = threshold))
}
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

CollectGarbage = function(){
  #if (exists("collect.garbage") ) rm(collect.garbage)
  ## The following function collects garbage until the memory is clean.
  ## Usage: 1. immediately call this function after you call a function or
  ##        2. rm()
  while (gc()[2,4] != gc()[2,4]){}
}

## creat cross-validation splits: this code is from the caret package 
createFolds <- function(y, k = 10, list = TRUE, returnTrain = FALSE) {  
  if(is.numeric(y)) {
    ## Group the numeric data based on their magnitudes
    ## and sample within those groups.
    
    ## When the number of samples is low, we may have
    ## issues further slicing the numeric data into
    ## groups. The number of groups will depend on the
    ## ratio of the number of folds to the sample size.
    ## At most, we will use quantiles. If the sample
    ## is too small, we just do regular unstratified
    ## CV
    cuts <- floor(length(y)/k)
    if(cuts < 2) cuts <- 2
    if(cuts > 5) cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  
  if(k < length(y)) {
    ## reset levels so that the possible levels and 
    ## the levels in the vector are the same
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    
    ## For each class, balance the fold allocation as far 
    ## as possible, then resample the remainder.
    ## The final assignment of folds is also randomized. 
    for(i in 1:length(numInClass)) {
      ## create a vector of integers from 1:k as many times as possible without 
      ## going over the number of samples in the class. Note that if the number 
      ## of samples in a class is less than k, nothing is producd here.
      min_reps <- numInClass[i] %/% k
      if(min_reps > 0) {
        spares <- numInClass[i] %% k
        seqVector <- rep(1:k, min_reps)
        ## add enough random integers to get  length(seqVector) == numInClass[i]
        if(spares > 0) seqVector <- c(seqVector, sample(1:k, spares))
        ## shuffle the integers for fold assignment and assign to this classes's data
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      } else {
        ## Here there are less records in the class than unique folds so
        ## randomly sprinkle them into folds. 
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k, size = numInClass[i])
      }  
    }
  } else foldVector <- seq(along = y)
  
  if(list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
    if(returnTrain) out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  } else out <- foldVector
  out
}

prettySeq <- function (x) paste("Resample", gsub(" ", "0", format(seq(along = x))), sep = "")

createResample <- function(y, times = 10, list = TRUE)
{
  trainIndex <- matrix(0, ncol = times, nrow = length(y))   
  out <- apply(trainIndex, 2, 
               function(data)
               {    
                 index <- seq(along = data)
                 out <- sort(sample(index, size = length(index), replace = TRUE))
                 out      
               })

  if (list) 
    {
      out <- as.data.frame(out)
      attributes(out) <- NULL
      names(out) <- prettySeq(out)
    } else {
      colnames(out) <- prettySeq(1:ncol(out))
    }
  
  out
}


cvSplit <- function(Nt, K){
size = Nt%/%K
A = runif(Nt)
rk = rank(A)
bk = factor((rk-1)%/%size + 1)
if(Nt%%K > 0) {
levels(bk)[levels(bk) == K+1]  = 1 
}
return(bk)
}




prediction_strength = function(data, min_num_clusters = 1, max_num_clusters = 10, 
                               num_trials = 5, ntree=20, no.rep = 10, thresh = 0.8, ...) {
	num_clusters = min_num_clusters:max_num_clusters
srt <- proc.time() 	
	prediction_strengths = c()
	for (i in 1:num_trials) {
		y = maply(num_clusters, function(n) 
		calculate_prediction_strength(data=data, num_clusters = n, ntree = ntree, no.rep = no.rep, ...))
		prediction_strengths = cbind(prediction_strengths, y)
	}
	
	means = aaply(prediction_strengths, 1, mean)
	stddevs = aaply(prediction_strengths, 1, sd)
	
# print(plot_prediction_strength(means, stddevs, num_clusters))
# We use 0.8 as our prediction strength threshold. Find the largest number of clusters with a 
# prediction strength greater than this threshold; this forms our estimate of the number of clusters.
	if (any(means > thresh)) {
		estimated = max((1:length(means))[means > thresh])		
	} else {
		estimated = which.max(means)
	}
ed <- proc.time()
run.time <- as.numeric((ed - srt)[3]) 
cat("Estimated run time = ", run.time, "\n") 
cat("====================================================== \n")	
cat("The estimated number of clusters =  ", estimated, "\n")	
return(estimated)
}



calculate_prediction_strength = function(data, num_clusters, ntree, no.rep, ...) {
	if (num_clusters == 1) {
		1 # The prediction strength is always 1 for 1 cluster.
	} else {
ix.cv <- cvSplit(Nt = nrow(data), K = 2)
training_set <- data[ix.cv == 1, ]
test_set <- data[ix.cv == 2, ] 
#		rands = runif(nrow(data), min = 0, max = 1)
#		training_set = data[(1:length(rands))[rands <= 0.5], ]
#		test_set = data[(1:length(rands))[rands > 0.5], ]		
		
d.trn <- RFdist(data= training_set, ntree = ntree, no.rep= no.rep, syn.type = "permute", importance= TRUE, ...)
d.tst <- RFdist(data= test_set, ntree = ntree, no.rep= no.rep, syn.type = "permute", importance= TRUE, ...)
hclust.trn <- hclust(d=d.trn$RFdist, method="ward.D2", members= NULL)
hclust.tst <- hclust(d=d.tst$RFdist, method="ward.D2", members= NULL)

trn.clusters <- cutree(hclust.trn, num_clusters)
tst.clusters <- cutree(hclust.tst, num_clusters)    
trn.centers <- mediod(x = d.trn$RFdist, clusters=trn.clusters)
trn.centers$mediod <- training_set[trn.centers$ids, ,drop = FALSE]		
# The prediction strength is the minimum prediction strength among all clusters.
prediction_strengths = maply(1:num_clusters, function(n) 
                       prediction_strength_of_cluster(test_set, tst.clusters, trn.centers, n))
min(prediction_strengths)
	}	
}


# Calculate the proportion of pairs of points in test cluster `k` that would again be 
# assigned to the same cluster, if each were clustered according to its closest training cluster mean.
prediction_strength_of_cluster = function(test_set, tst.clusters, trn.centers, k) {
	if (sum(tst.clusters == k) <= 1) {
		1 # No points in the cluster.
	} else {
		test_cluster = test_set[tst.clusters == k, ,drop = FALSE]
		count = 0
		for (i in 1:(nrow(test_cluster)-1)) {
			for (j in (i+1):nrow(test_cluster)) {
				p1 = data.matrix(test_cluster[i, , drop = FALSE])
				p2 = data.matrix(test_cluster[j, , drop = FALSE])			
                        l1 <- predict(object=trn.centers, org.data=NULL, newdata = p1)
                        l2 <- predict(object=trn.centers, org.data=NULL, newdata = p2)
	                if (l1 == l2) {
				count = count + 1
				}
			}
		}
		# Return the proportion of pairs that stayed in the same cluster.
		count / (nrow(test_cluster) * (nrow(test_cluster) - 1) / 2) 
	}
}

plot_prediction_strength = function(means, stddevs, num_clusters) {
	qplot(num_clusters, means, xlab = "# clusters", ylab = "prediction strength", geom = "line", 
	main = "Estimating the number of clusters via the prediction strength") + 
	geom_errorbar(aes(num_clusters, ymin = means - stddevs, ymax = means + stddevs), 
	size = 0.3, width = 0.2, colour = "darkblue")
}

multiClassSummary <- function (data, lev = NULL, model = NULL)
{
 #   levels(data[, "pred"]) = lev 
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
        stop("levels of observed and predicted data do not match")
        
    has_class_probs <- all(lev %in% colnames(data))
    if (has_class_probs) {
        lloss <- mnLogLoss(data = data, lev = lev, model = model)
        requireNamespaceQuietStop("ModelMetrics")
        prob_stats <- lapply(levels(data[, "pred"]), function(x) {
            obs <- ifelse(data[, "obs"] == x, 1, 0)
            prob <- data[, x]
            AUCs <- try(ModelMetrics::auc(obs, data[, x]), silent = TRUE)
            return(AUCs)
        })
        roc_stats <- mean(unlist(prob_stats))
    }
    CM <- confusionMatrix(data[, "pred"], data[, "obs"])
    if (length(levels(data[, "pred"])) == 2) {
        class_stats <- CM$byClass
    }
    else {
        class_stats <- colMeans(CM$byClass)
        names(class_stats) <- paste("Mean", names(class_stats))
    }
    overall_stats <- if (has_class_probs)
        c(CM$overall, logLoss = lloss, ROC = roc_stats)
    else CM$overall
    if (length(levels(data[, "pred"])) > 2)
        names(overall_stats)[names(overall_stats) == "ROC"] <- "Mean_AUC"
    stats <- c(overall_stats, class_stats)
    stats <- stats[!names(stats) %in% c("AccuracyNull", "AccuracyLower",
        "AccuracyUpper", "AccuracyPValue", "McnemarPValue", "Mean Prevalence",
        "Mean Detection Prevalence")]
    names(stats) <- gsub("[[:blank:]]+", "_", names(stats))
    stat_list <- c("Accuracy", "Kappa", "Mean_Sensitivity", "Mean_Specificity",
        "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", "Mean_Detection_Rate",
        "Mean_Balanced_Accuracy")
    if (has_class_probs)
        stat_list <- c("logLoss", "Mean_AUC", stat_list)
    if (length(levels(data[, "pred"])) == 2)
        stat_list <- gsub("^Mean_", "", stat_list)
    stats <- stats[c(stat_list)]
    return(stats)
}


comembership_table <- function(labels1, labels2) {
  if (length(labels1) != length(labels2)) 
    stop("The two vectors of cluster labels must be of equal length.");
    .Call("rcpp_comembership_table", labels1, labels2, PACKAGE = "UnsupRF")
}


calinhara <- function (x, clustering, cn = max(clustering)) {
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    cln <- rep(0, cn)
    W <- matrix(0, p, p)
    for (i in 1:cn) cln[i] <- sum(clustering == i)
    for (i in 1:cn) {
        clx <- x[clustering == i, ]
        cclx <- cov(as.matrix(clx))
        if (cln[i] < 2)
            cclx <- 0
        W <- W + ((cln[i] - 1) * cclx)
    }
    S <- (n - 1) * cov(x)
    B <- S - W
    a = (cn - 1)/(n - cn)
    SSB <-  sum(diag(B))
    (a/SSB)*diag(W)
}

##
Cluster <- function(data, k, cluster.method, control){
#arg <- list(...) 
#nme <- names(arg) 
  switch(cluster.method, 
         UnsupRF={
           distance <- RFdist(data=data,  syn.type = "permute", importance= FALSE, no.rep = control$no.rep, 
           ntree = control$ntree, parallel = control$RF.parallel)$RFdist
           clust.model <- hclust(d=distance, method=control$method, members= NULL)
           clusters <- cutree(clust.model, k)   
         }, 
         kmeans = {
         distance <- dist(data) 
#           if( ("iter.max"%in%nme) | ("algorithm"%in%nme) |  ("nstart" %in% nme)) {
#           iter.max = arg$iter.max; algorithm = arg$algorithm;  nstart = arg$nstart 
#           clusters <- kmeans(x = data, centers = k, iter.max= iter.max, algorithm = algorithm, nstart = nstart)$cluster
           clusters <- kmeans(x = data, centers = k)$cluster
         }, 
         hclust = {
           distance <- dist(data)
           clust.model <- hclust(d=distance, method=control$method, members= NULL)
           clusters <- cutree(clust.model, k)   
         }, 
         stop("Wrong value for cluster method")
  )
  list(clusters = clusters, distance = distance) 
}













