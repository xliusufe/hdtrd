pvaltrans_cv <- function(target, source, family = "gaussian", delta0 = 0.1, nsource = 10, method = 'lasso', ncv = 10, alpha = 0.05, resids = NULL, isproj = FALSE, proj = NULL, sigma2 = NULL, lammax = NULL, nmoms = NULL, zK = NULL, J = NULL, K = NULL, timeout = 0){
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

	X0	= target$X
	Y0 	= target$Y
	n 	= nrow(X0)
	p 	= ncol(X0)

	lammax_all 	= rep(0, nsource)
	for(k in 1:nsource){
		datsource = source[[k]]
		X 	= datsource$X
		ns 	= nrow(X)
		if(ns>n){
			Xk = X[1:n,]
			Xd = X0-Xk
		}
		else{
			Xk = X
			Xd = X0[1:ns,]-Xk
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
			Xp 	= as.matrix(projk)
		}
		else{
			Xp = Xk
		}

		K	= ifelse(is.null(K), max(500,2*n,2*p), K)
		J	= ifelse(is.null(J), max(500,2*n,2*p), J)
		if(lammax=='mpmo'){
			lammax1 <- eigmax(X = Xp, zK = zK, J = J, K = K, method = 'mpmo', nmoms = nmoms, timeout = timeout)$lammax
		}
		else if(lammax=='mplp'){
			lammax1 <- eigmax(X = Xp, zK = zK, J = J, K = K, method = 'mplp', nmoms = nmoms, timeout = timeout)$lammax
		}
		else if(lammax=='empi'){
			lammax1 <- eigmax(X = Xp, zK = zK, J = J, K = K, method = 'empi', nmoms = nmoms, timeout = timeout)$lammax
		}
		else{
			stop("lammax must be a numeric number or be on of {'mpmp','mplp','empi'} !")
		}
		lammax_all[k] = lammax1
	}


	ndelta 	= length(delta0)
	RSS 	= rep(0, ndelta)
	len_cv	= ceiling(n/ncv)
	for(jj in 1:ncv){
		cv.id	= ((jj-1)*len_cv+1):(jj*len_cv)
		if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
		ytrain	= Y0[-cv.id]
		xtrain	= X0[-cv.id,]
		ytest	= Y0[cv.id]
		xtest	= X0[cv.id,]
		nt 		= length(ytrain)

		pvals	= matrix(0, nsource, ndelta)
		for(k in 1:nsource){
			datsource = source[[k]]
			X 	= datsource$X
			Y 	= datsource$Y

			ns 	= length(Y)
			if(ns>nt){
				Xk = X[1:nt,]
				Xd = xtrain-Xk
				Yd = ytrain-Y[1:nt]
			}
			else{
				Xk = X
				Xd = xtrain[1:ns,]-Xk
				Yd = ytrain[1:ns]-Y
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
				Xp 	= as.matrix(projk[-cv.id,])
			}
			else{
				Xp = Xk[-cv.id,]
			}


			if(is.numeric(lammax)||is.null(lammax)){
				Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0, method = method, resids = resids, sigma2 = sigma2, lammax = lammax)
			}
			else{
				Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0, method = method, resids = resids, sigma2 = sigma2, lammax = lammax_all[k])
			}
			pvals[k,]	= Tn1$pvals
		}

		for(kk in 1:ndelta){
			idtrans 	<- which(pvals[, kk]>alpha)
			fittrans	<- utrans(target=list(X=xtrain,Y=ytrain), source=source, family = family, idtrans=idtrans)
			mu			<- fittrans$beta[1] + xtest%*%fittrans$beta[-1]
			if(family == 'gaussian'){
				pred_err = mean((ytest - mu)^2)
			}
			else if(family == 'binomial'){
				pred_err = mean((ytest - 1.0/(1.0+exp(-mu)))^2)
			}
			else if(family == 'poisson'){
				pred_err = mean((ytest - exp(mu))^2)
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}
			RSS[kk] = RSS[kk] + pred_err
		}
	}
	s_opt = which.min(RSS)

	pvals	= rep(0, nsource)
	for(k in 1:nsource){
		datsource = source[[k]]
		X 	= datsource$X
		Y 	= datsource$Y

		ns 	= length(Y)
		if(ns>n){
			Xk = X[1:n,]
			Xd = X0-Xk
			Yd = Y0-Y[1:n]
		}
		else{
			Xk = X
			Xd = X0[1:ns,]-Xk
			Yd = Y0[1:ns]-Y
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
			Xp 	= as.matrix(projk)
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

		if(is.numeric(lammax)||is.null(lammax)){
			Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0[s_opt], method = method, resids = resids, sigma2 = sigma2, lammax = lammax)
		}
		else{
			Tn1	= pvalrd(data = list(Y=Yd,X=Xp,Z=Xd), family = family, delta0 = delta0[s_opt], method = method, resids = resids, sigma2 = sigma2, lammax = lammax_all[k] )
		}
		pvals[k]	= Tn1$pvals
	}

	return(list(pvals=pvals, s_opt=s_opt))
}
