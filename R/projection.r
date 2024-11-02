projection <- function(x, y, family = "gaussian", method = 'lasso', isresid = TRUE){
	p 	= ifelse(is.null(ncol(y)), 1, ncol(y))
	n 	= ifelse(is.null(nrow(x)), 1, length(x))
	proj = NULL
	if(method == "lasso"){
		if(is.null(ncol(y))){
			fit_k	= cv.glmnet(x, y, type.measure="mse")
			proj	= cbind(proj, coef(fit_k))
		}
		else{
			for(k in 1:p){
				fit_k	= cv.glmnet(x, y[,k], type.measure="mse")
				proj	= cbind(proj, coef(fit_k))
			}
		}
	}
	else if(method == "glm"){
		p1 	= ifelse(is.null(ncol(x)), 1, ncol(x))
		if(p1>n/2){
			stop("Dimension of X is greater than n/2! Please use mathod = lasso.")
		}
		if(is.null(ncol(y))){
			fit_k 	= glm(y~x, family = family)
			proj 	= cbind(proj, coef(fit_k))
		}
		else{
			for(k in 1:p){
				fit_k 	= glm(y[,k]~x, family = family)
				proj 	= cbind(proj, coef(fit_k))
			}
		}
	}
	else{
		stop("Method must be one of {'glm', 'lasso'} !")
	}
	if(isresid){
		proj = y - cbind(1,x)%*%proj
	}

	return(proj)
}