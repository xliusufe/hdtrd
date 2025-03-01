translasso0 <- function(target, source = NULL, idtrans = NULL, nvec = NULL, lamconst = NULL, l1 = TRUE){
	if(is.null(target)){
		stop("target dataset can not be NULL !")
	}
	if(is.null(source)){
		X 	= target$X
		Y 	= target$Y
		p	<- ncol(X)
		if(is.null(nvec)){
			nvec = length(target$Y)
		}
		cv.init	<- cv.glmnet(X[1:nvec[1],], Y[1:nvec[1]], nfolds=8, lambda=seq(1,0.1,length.out=20)*sqrt(2*log(p)/nvec[1]))
		lamconst<- cv.init$lambda.min/sqrt(2*log(p)/nvec[1])
		beta.kA <- predict(cv.init, s='lambda.min', type='coefficients')[-1]
		w.kA<-NA
		return(list(beta	= c(0, as.numeric(beta.kA)),
					w.kA	= w.kA)
			)
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
	p	<- ncol(X)
	size.idtrans <- length(idtrans)
	ind.kA	<- indset(nvec, c(1, idtrans+1))
	ind.1	<- 1:nvec[1]
	if(l1){
		Y.A	<- Y[ind.kA]
	}else{ #the l0-method
		Y.A<- Y[ind.1]
		Sig.hat<-t(X)%*%X/nrow(X)
		for(k in 1:size.idtrans){
			ind.k	<- indset(nvec,k+1)
			lam.k 	<- sqrt(mean(Y[ind.1]^2)/nvec[1]+mean(Y[ind.k]^2)/nvec[k]) * sqrt(2*log(p))
			delta.hat.k	<- glmnet(XtX	= Sig.hat, 
										Xty		= t(X[ind.k,])%*%Y[ind.k]/nvec[k+1]-t(X[1:nvec[1],])%*%Y[1:nvec[1]]/nvec[1],
										lambda	= lam.k)$coef
			Y.A	<- c(Y.A, Y[ind.k]-X[ind.k,]%*%delta.hat.k)
		}
	}

	if(is.null(lamconst)){
		cv.init	<- cv.glmnet(X[ind.kA,], Y.A, nfolds=8, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.kA)))
		lamconst<- cv.init$lambda.min/sqrt(2*log(p)/length(ind.kA))
	}
	w.kA	<- as.numeric(glmnet(X[ind.kA,], Y.A, lambda=lamconst*sqrt(2*log(p)/length(ind.kA)))$beta)
	w.kA	<- w.kA*(abs(w.kA)>=lamconst*sqrt(2*log(p)/length(ind.kA)))

	delta.kA	<- as.numeric(glmnet(x=X[ind.1,],y=Y[ind.1]-X[ind.1,]%*%w.kA, lambda=lamconst*sqrt(2*log(p)/length(ind.1)))$beta)
	delta.kA	<- delta.kA*(abs(delta.kA)>=lamconst*sqrt(2*log(p)/length(ind.1)))
	beta.kA 	<- w.kA + delta.kA
	lamconst	= NA

	return(list(beta	= c(0, as.numeric(beta.kA)),
				w.kA	= w.kA)
		)
}