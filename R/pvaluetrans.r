pvaltrans <- function(target, source, family = "gaussian", delta0 = 0.1, nsource = 10, testmethd = 'pvalrd', method = 'lasso', resids = NULL, isproj = FALSE, proj = NULL, sigma2 = NULL, lammax = NULL, nmoms = NULL, zK = NULL, J = NULL, K = NULL, timeout = 0){
	# High-dimensional testing of relevant difference in linear regressions.
	if(is.null(target)){
		stop("target  dataset can not be NULL !")
	}
	if(is.null(source)){
		stop("sourcr  dataset can not be NULL !")
	}
	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	if(nsource<1){
		stop("the number of source datasets is not smaller than one!")
	}

	if(is.null(lammax)){
		lammax = "empi"
	}
	x	= target$X
	y 	= target$Y
	n 	= length(y)
	p 	= ncol(x)

	pvals	= rep(0, nsource)
	for(k in 1:nsource){
		datsource = source[[k]]
		X 	= datsource$X
		Y 	= datsource$Y
		
		if(testmethd == 'pvalrd'){
			ns 	= length(Y)
			if(ns>n){
				Xk = X[1:n,]
				Xd = x-Xk
				Yd = y-Y[1:n]
			}
			else{
				Xk = X
				Xd = x[1:ns,]-Xk
				Yd = y[1:ns]-Y
			}
			if(isproj){
				if(is.null(proj)){
					projk = projection(Xd, Xk, family = family, method = method)
				}
				else{
					if(is.list(proj)){
						projk = proj[[k]]
					}
					else if(is.array(proj)){
						projk = proj[,,k]
					}
				}
				# Xp 	= Xk - cbind(1,Xd)%*%projk
				Xp 	= projk
			}
			else{
				Xp = Xk
			}

			if(method == "lasso"){
				fitglm 	= cv.glmnet(Xd, Yd, family = family, type.measure="mse")
				betahat = coef(fitglm)
			}
			else if(method == "glm"){
				x 	= target$X
				p1 	= ifelse(is.null(ncol(x)), 1, ncol(x))
				if(p1>n/2){
					stop("Dimension of X is greater than n/2! Please use mathod = lasso.")
				}
				fitglm 	= glm(Yd~Xd, family = family)
				betahat = coef(fitglm)
			}
			else{
				stop("If target$X is null, method must be one of {'glm', 'lasso'} !")
			}

			mu 		= betahat[1] + Xd %*% betahat[-1]
			if(family=='gaussian'){
				resids	<- Yd - mu
			}
			else if(family == 'binomial'){
				resids	<- Yd - 1/(1+exp(-mu))
			}
			else if(family == 'poisson'){
				resids	<- Yd - exp(mu)
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}

		
			K	= ifelse(is.null(K), max(500,2*n,2*p), K)
			J	= ifelse(is.null(J), max(500,2*n,2*p), J)
			if(lammax=='mpmo'){
				lammax1 <- eigmax(X = Xp, zK = zK, J = J, K = K, method = 'mpmo', nmoms = nmoms, timeout = timeout)$lammax
				Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0, method = method, resids = resids, sigma2 = sigma2, lammax = lammax1)
			}
			else if(lammax=='mplp'){
				lammax1 <- eigmax(X = Xp, zK = zK, J = J, K = K, method = 'mplp', nmoms = nmoms, timeout = timeout)$lammax
				Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0, method = method, resids = resids, sigma2 = sigma2, lammax = lammax1)
			}
			else if(lammax=='empi'){
				lammax1 <- eigmax(X = Xp, zK = zK, J = J, K = K, method = 'empi', nmoms = nmoms, timeout = timeout)$lammax
				Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0, method = method, resids = resids, sigma2 = sigma2, lammax = lammax1)
			}
			else if(is.numeric(lammax)||is.null(lammax)){
				Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0, method = method, resids = resids, sigma2 = sigma2, lammax = lammax)
			}
			else{
				stop("lammax must be a numeric number or be on of {'mpmp','mplp','empi'} !")
			}
		}
		else{
			X01 = rbind(x,X)
			Y01 = c(y,Y)
			fitglm 	= cv.glmnet(X01, Y01, family = family, type.measure="mse")
			betahat = coef(fitglm)
			mu 		= betahat[1] + X01 %*% betahat[-1]
			resids	<- Y01 - mu
			Tn1		= pvalclc(data = list(Y=Y01, X=X01,Z=X), family = family, method = method, resids = resids)
			# Tn1	= pvalclc(data = list(Y=Yd,X=Xp,Z=Xd), family = family, method = method, resids = resids)
		}
		pvals[k]	= Tn1$pvals
	}

	return(pvals)
}
