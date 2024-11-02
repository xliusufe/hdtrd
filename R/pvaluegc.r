pvalgc <- function(data, family = "gaussian", resids = NULL, psi = NULL){
	# High-dimensional testing of coefficient in linear regressions.
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
		else{
			z 		= data$Z
			fitglm 	= glm(y~z, family = family)
			resids 	= residuals(fitglm, type = "response")
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
	pvals 	= pnorm(Tn, lower.tail = F)

	return(list(Tn = Tn, pvals = pvals))
}
