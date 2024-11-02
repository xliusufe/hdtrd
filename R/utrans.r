utrans <- function(target, source = NULL, family = "gaussian", idtrans = NULL){
	if(is.null(target)){
		stop("target dataset can not be NULL !")
	}
	if(is.null(source)){
		fitglm 	= cv.glmnet(target$X, target$Y, family = family, type.measure="mse")
		betahat = coef(fitglm)
		return(list(fitglmnet = fitglm, beta = betahat, family = family))
	}

	nsource = length(source)
	if(is.null(idtrans)){
		idtrans = seq(nsource)
	}	
	else if(length(setdiff(idtrans,seq(nsource)))>0){
		stop("idtrans must be a subset of 1,...,ns!")
	}

	nid = length(idtrans)

	y 	= target$Y
	x0 	= target$X
	p 	= ncol(x0)

	if(nid<=0){
		x = x0
	}
	else{
		xs = list()
		xs[[1]] = x0
		for(k in 1:nid){
			dsource_k = source[[idtrans[k]]]
			y	= c(y, dsource_k$Y )
			x0 	= rbind(x0, dsource_k$X)
			xs[[k+1]] = dsource_k$X
		}
		x = bdiag(xs)
		x[,1:p] = x0
	}

	fitglm 	= cv.glmnet(x, y, family = family, type.measure="mse")
	betahat = coef(fitglm)[1:(p+1)]
	return(list(fitglmnet = fitglm, beta = betahat, family = family))
}