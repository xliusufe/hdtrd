pvalrd <- function(data, family = "gaussian", delta0 = 0.1, method = 'lasso', resids = NULL, sigma2 = NULL, lammax = NULL){
	# High-dimensional testing of relevant difference in linear regressions.

	multidelta = 0
	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	if(is.null(data$X)){
		stop("Z is NULL !")
	}
	if(is.null(data$Y)){
		stop("Y is NULL !")
	}
	y 	= data$Y
	n 	= length(y)
	x 	= data$X
	p 	= ifelse(is.null(ncol(z)), 1, ncol(z))

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
			# resids 	= y - predict(fitglm, s = "lambda.min", newx = x)
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
		else if(method == "glm"){
			z 	= data$Z
			p1 	= ifelse(is.null(ncol(z)), 1, ncol(z))
			if(p1>n/2){
				stop("Dimension of X is greater than n/2! Please use mathod = lasso.")
			}
			fitglm 	= glm(y~z, family = family)
			resids 	= residuals(fitglm, type = "response")
		}
		else{
			stop("If data$X is null, method must be one of {'glm', 'lasso'} !")
		}
	}

	ishd = 1
	if(is.null(lammax)){
		lammax = norm(t(x)%*%x/n,"2")/(1+sqrt(p/n))
	}
	if(ishd){
		maxsig = ( delta0*lammax )^2
	}
	else{
		maxsig = delta0^2*lammax
	}

	if(is.null(sigma2)){
		sigma2 	= 1
	}

	if(multidelta == 1){
		dims 	= c(n, p, length(delta0), ishd)
		Tn		= .Call("RDtest_MAX_GC_CV1_",
						as.numeric(x),
						as.numeric(resids),
						as.numeric(sigma2),
						as.numeric(maxsig),
						as.integer(dims)
				)
	}
	if(multidelta == 2){
		dims 	= c(n, p, length(delta0), ncol(resids), ishd)
		Tn		= .Call("RDtest_MAX_GC_CV2_",
						as.numeric(x),
						as.numeric(resids),
						as.numeric(sigma2),
						as.numeric(maxsig),
						as.integer(dims)
				)
	}
	else{
		dims 	= c(n, p, ishd)
		Tn		= .Call("RDtest_MAX_GC_",
						as.numeric(x),
						as.numeric(resids),
						as.numeric(sigma2),
						as.numeric(maxsig),
						as.integer(dims)
				)
	}
	pvals 	= pnorm(Tn, lower.tail = F)

	return(list(Tn = Tn, pvals = pvals))
}
