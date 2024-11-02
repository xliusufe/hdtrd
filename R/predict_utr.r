predict_utr <- function(fittrans, X, type = "response"){
	familyk = fittrans$family
	mu = fittrans$beta[1] + X%*%fittrans$beta[-1]
	if(familyk == 'gaussian'){
		yhat = mu
	}
	else if(familyk == 'binomial'){
		if(type == "class"){
			yhat = as.numeric(mu>=0)
		}
		else if(type == "response"){
			yhat = 1.0/(1.0 + exp(-mu))
		}
	}
	else if(familyk == 'poisson'){
		yhat = exp(mu)
	}
	else{
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	return(yhat)
}