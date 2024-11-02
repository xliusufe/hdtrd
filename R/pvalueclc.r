pvalclc <- function(data, family = "gaussian", method = 'lasso', resids = NULL, psi = NULL){
	# High-dimensional testing of coefficient in linear regressions in presence of high-dimensional control factors.
	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	if(is.null(data$X)){
		stop("X is NULL !")
	}
	if(is.null(data$Y)){
		stop("Y is NULL !")
	}
	y 	= data$Y
	n 	= length(y)
	x 	= data$X
	p 	= ifelse(is.null(ncol(x)), 1, ncol(x))

	if(is.null(resids)){
		if(is.null(data$Z)){
			if(family=='gaussian'){
				resids  = y
			}
			else if(family == 'binomial'){
				resids  = y - 0.5
			}
			else if(family == 'poisson'){
				resids  = y - 1
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}
		}
		else if(method == "lasso"){
			z 		= data$Z
			fitglm 	= cv.glmnet(z, y, family = family, type.measure="mse")
			betahat = coef(fitglm)
			mu 		= betahat[1] + z %*% betahat[-1]

			if(family=='gaussian'){
				resids	<- y - mu
			}
			else if(family == 'binomial'){
				resids	<- y - 1/(1+exp(-mu))
			}
			else if(family == 'poisson'){
				resids	<- y - exp(mu)
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}
		}
		else{
			z 		= data$Z
			fitglm	<- gfabs(list(Y = y, X = z), family="gaussian")
			resids	= fitglm$residual
		}
	}
	if(is.null(psi)){
		psi     = rep(1, n)
		ispsi   = 0
	}
	else{
		ispsi = 1
	}


	dims 	= c(n, p, ispsi)
	Tn		= .Call("GCtest_",
				as.numeric(x),
				as.numeric(resids),
				as.numeric(psi),
				as.integer(dims)
			)
	pvals 	= pnorm(abs(Tn), lower.tail = F)

	return(list(Tn = abs(Tn), pvals = pvals))
}
