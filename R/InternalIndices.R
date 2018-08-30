

#' @export
index.DB<-function(x,cl,d=NULL,centrotypes="centroids",p=2,q=2){
   	 if(sum(c("centroids","medoids")==centrotypes)==0)
      stop("Wrong centrotypes argument")
   	 if("medoids"==centrotypes && is.null(d))
      stop("For argument centrotypes = 'medoids' d cannot be null")
     if(!is.null(d)){
      if(!is.matrix(d)){
        d<-as.matrix(d)
      }
     row.names(d)<-row.names(x)
     }
   if(is.null(dim(x))){
          dim(x)<-c(length(x),1)
        }
  x<-as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  #print(n)
  dAm<-d
  centers<-matrix(nrow=k,ncol=ncol(x))
  if (centrotypes=="centroids"){
    for(i in 1:k)
    {
      for(j in 1:ncol(x))
      {
         centers[i,j]<-mean(x[cl==i,j])
      }
    }
  }
  else if (centrotypes=="medoids"){
    #print("start")
    #print(dAm)
    for (i in 1:k){
      clAi<-dAm[cl==i,cl==i]
      if (is.null(clAi)){
        centers[i,]<-NULL
      }
      else{
        #print("przed centers")
        #print(x[cl==i,])
        #print(clAi)
        centers[i,]<-.medoid(x[cl==i,],dAm[cl==i,cl==i])
        #print("po centers")
        #print(centers[i])
      }
    }   
    #print("stop")
  }
  else{
    stop("wrong centrotypes argument")
  }
  S<-rep(0,k)
  for(i in 1:k){                             # For every cluster
    ind <- (cl==i)
    if (sum(ind)>1){
      centerI<-centers[i,]
      centerI<-rep(centerI,sum(ind))
      centerI<-matrix(centerI,nrow=sum(ind),ncol=ncol(x),byrow=TRUE)
      S[i] <- mean(sqrt(apply((x[ind,] - centerI)^2,1,sum))^q)^(1/q)
    }
    else
      S[i] <- 0                         
  }
  M<-as.matrix(dist(centers,p=p))
  R <- array(Inf,c(k,k))
  r = rep(0,k)
  for (i in 1:k){
    for (j in 1:k){
        R[i,j] = (S[i] + S[j])/M[i,j]
    }
    r[i] = max(R[i,][is.finite(R[i,])])
  } 
  DB = mean(r[is.finite(r)])        
  resul<-list(DB=DB,r=r,R=R,d=M,S=S,centers=centers)
  resul
}

#' @exportindex.G2 <- function(d,cl){
  cn <- max(cl)
  n <- length(cl)
  dmat <- as.matrix(d)
  diameter <- average.distance <- median.distance <- separation <-
    average.toother <- 
    cluster.size <- within.dist <- between.dist <- numeric(0)
  separation.matrix <- matrix(0,ncol=cn,nrow=cn)
  di <- list()
  for (i in 1:cn){
    cluster.size[i] <- sum(cl==i)
    #print(i)
    #print(cl==i)
    #print(dmat[cl==i,cl==i])
    di <- as.dist(dmat[cl==i,cl==i])
    within.dist <- c(within.dist,di)
    #diameter[i] <- max(di)
    average.distance[i] <- mean(di)
    median.distance[i] <- median(di)
    bv <- numeric(0)
    for (j in 1:cn){
      if (j!=i){
        sij <- dmat[cl==i,cl==j]
        bv <- c(bv,sij)
        if (i<j){
          separation.matrix[i,j] <- separation.matrix[j,i] <- min(sij)
          between.dist <- c(between.dist,sij)
        }
      }
    }   
   }
   nwithin<-length(within.dist)
   nbetween<-length(between.dist)
   .C("fng2",as.double(within.dist),as.integer(nwithin),as.double(between.dist),as.integer(nbetween),wynik=double(1),PACKAGE="UnsupRF")$wynik[1]
}

#' @export
index.G3<-function(d,cl)
{
	d<-data.matrix(d)
     	DU<-0
	r<-0
	v_max<-array(1,max(cl))
	v_min<-array(1,max(cl))
	for (i in 1:max(cl))
	{
		
		n<-sum(cl==i)
		if (n>1)
		{
			t<-d[cl==i,cl==i]
			DU=DU+sum(t)/2
			v_max[i]=max(t)
			if (sum(t==0)==n)		# nie ma zer poza przekatna
				v_min[i]<-min(t[t!=0])
			else
				v_min[i]<-0
			r<-r+n*(n-1)/2
		}
	}
	Dmin=min(v_min)
	Dmax=max(v_max)
	if(Dmin==Dmax)
		result<-NA
	else
		result<-(DU-r*Dmin)/(Dmax*r-Dmin*r)
	result		
}


