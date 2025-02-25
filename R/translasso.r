indset <- function(nvec, k.vec){
	ind.re <- NULL
	for(k in k.vec){
		if(k==1){
			ind.re<-c(ind.re,1: nvec[1])
		}else{
			ind.re<- c(ind.re, (sum(nvec[1:(k-1)])+1): sum(nvec[1:k]))
		}
	}
	return(ind.re)
}

agg.fun <- function(B, X.test,y.test, total.step=10, selection=F){
	if(sum(B==0)==ncol(B)*nrow(B)){
		return(rep(0,nrow(B)))
	}
	p<-nrow(B)
	K<-ncol(B)
	colnames(B)<-NULL
	if(selection){#select beta.hat with smallest prediction error
		khat<-which.min(colSums((y.test-X.test%*%B)^2))
		theta.hat<-rep(0, ncol(B))
		theta.hat[khat] <- 1
		beta=B[,khat]
		beta.ew=NULL
	}else{#Q-aggregation
		theta.hat	<- exp(-colSums((y.test-X.test%*%B)^2)/2)
		theta.hat	= theta.hat/sum(theta.hat)
		theta.old	= theta.hat
		beta		<- as.numeric(B%*%theta.hat)
		beta.ew		<- beta
		for(ss in 1:total.step){
			theta.hat	<- exp(-colSums((y.test-X.test%*%B)^2)/2+colSums((as.vector(X.test%*%beta)-X.test%*%B)^2)/8)
			theta.hat	<- theta.hat/sum(theta.hat)
			beta		<- as.numeric(B%*%theta.hat*1/4+3/4*beta)
			if(sum(abs(theta.hat-theta.old))<10^(-3)){
				break
			}
			theta.old=theta.hat
		}
	}
	return(list(theta=theta.hat, beta=beta, beta.ew=beta.ew))
}

translasso <- function(target, source = NULL, idtrans = NULL, nvec = NULL, Itil = NULL, l1 = TRUE){
	if(is.null(target)){
		stop("target dataset can not be NULL !")
	}
	if(is.null(source)){
		fittrans = translasso0(target, source = NULL)
		return(fittrans)
	}

	nsource = length(source)
	if(is.null(idtrans)){
		idtrans = seq(nsource)
	}	
	else if(length(setdiff(idtrans,seq(nsource)))>0){
		stop("idtrans must be a subset of 1,...,ns!")
	}

	X = target$X
	Y = target$Y
	for(k in 1:nsource)	{
		X	= rbind(X, source[[k]]$X)
		Y 	= c(Y, source[[k]]$Y)
	}

	if(is.null(nvec)){
		nvec = length(target$Y)
		for(k in 1:nsource)	{
			nvec = c(nvec, length(source[[k]]$Y))
		}
	}

	if(is.null(Itil)){
		n0 	= length(target$Y)
		Itil = c(1:min(50, ceiling(n0/3)))
	}

	M		= length(nvec)-1
	#step 1
	X0.til	<- X[Itil,] #used for aggregation
	y0.til	<- Y[Itil]
	X		<- X[-Itil,]
	y		<- Y[-Itil]

	#step 2
	Rhat	<- rep(0, M+1)
	p		<- ncol(X)
	nvec[1]	<- nvec[1]-length(Itil)
	ind.1	<- indset(nvec,1)
	for(k in 2: (M+1)){
		ind.k	<- indset(nvec,k)
		Xty.k 	<- t(X[ind.k,])%*%y[ind.k]/nvec[k] - t(X[ind.1,])%*%y[ind.1]/ nvec[1]
		margin.T<- sort(abs(Xty.k),decreasing=T)[1:round(nvec[1]/3)]
		Rhat[k] <- sum(margin.T^2)
	}
	Tset 	<- list()
	k0		=0
	kk.list	<-unique(rank(Rhat[-1]))
	for(kk in 1:length(kk.list)){#use Rhat as the selection rule
		Tset[[k0+kk]]	<- which(rank(Rhat[-1]) <= kk.list[kk])
	}
	k0	=length(Tset)
	Tset<- unique(Tset)

	target_Itil = list(X = target$X[-Itil,], Y = target$Y[-Itil])
	beta.T	<- list()
	init.re	<- translasso0(target=target_Itil, source=source, idtrans=NULL, nvec=nvec, l1=l1)
	beta.T[[1]] 	<- init.re$beta[-1]
	beta.pool.T		<- beta.T ##another method for comparison
	for(kk in 1:length(Tset)){#use pi.hat as selection rule
		re.k	<- translasso0(target=target_Itil, source=source, idtrans=Tset[[kk]], nvec=nvec, l1=l1, lamconst=init.re$lam.const)
		beta.T[[kk+1]] 		<-re.k$beta[-1]
		beta.pool.T[[kk+1]]	<-re.k$w.kA
	}
	beta.T	<- beta.T[!duplicated((beta.T))]
	beta.T	<- as.matrix(as.data.frame(beta.T))
	agg.re1	<- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til)
	beta.pool.T	<- beta.pool.T[!duplicated((beta.pool.T))]
	beta.pool.T	<- as.matrix(as.data.frame(beta.pool.T))
	agg.re2		<- agg.fun(B=beta.pool.T, X.test=X0.til, y.test=y0.til)

	return(list(beta		= c(0, as.numeric(agg.re1$beta)), 
				rankpi		= rank(Rhat[-1]))
	)
}