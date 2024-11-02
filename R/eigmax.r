eigmax <- function(X, zK = NULL, tJ = NULL, K = 1000, J = 1000, method = 'mpmo', nmoms = NULL, timeout = 0){

	n 	= nrow(X)
	p 	= ifelse(is.null(ncol(X)), 1, ncol(X))

	tau	= p/n
	Sn	= t(X)%*%X/n
	if(method=='mplp'){
		if(is.null(nmoms)){
			nmoms = 0
		}
		else if(nmoms>4 || nmoms <0 ){
			nmoms =  0
			warning("Mumber of moments is in {0, 1, 2, 3, 4}!")
		}

		if(is.null(zK)){
			zr = rnorm(K)
			zi = rep(1,K)/sqrt(n)
		}
		else{
			zr	= zK[,1] # real part
			zi	= zK[,2] # imaginary part
			kz 	= nrow(zK)
			if(K != kz){
				options(warn=1)
				warning("The number of columns of Z is not K. We reset K to be the number of columns of Z.")
				K = kz
			}
		}

		eigs	= eigen(Sn)$values
		fitmom	= .Call("_ARGUEMENTS_MPLP4",
					as.numeric(zr),
					as.numeric(zi),
					as.numeric(eigs),
					as.integer(c(n, p, nmoms, J, K))
					)

		tJ = fitmom$tJ
		moms = fitmom$moms

		fit	<- try(
				linp( E	= t(matrix(fitmom$E, ncol = nmoms+1)), 
					G 	= t(matrix(fitmom$G, nrow = J+2*K)), 
					F 	= fitmom$F,
					H 	= fitmom$H,
					Cost = fitmom$f,
					verbose = FALSE,
					timeout = timeout
				),
			silent=FALSE
			)

		if('try-error' %in% class(fit) || fit$IsError){
			Sn		= X%*%t(X)/n
			eigs	= eigen(Sn)$values

			fitmom	= .Call("_ARGUEMENTS_MPMOM",
						as.numeric(Sn),
						as.numeric(eigs),
						as.integer(c(n, p, 4, J, 0))
						)
			tJ = fitmom$tJ

			fit	<- try(
					linp( E 	= fitmom$E,
						F 	= 1,
						G 	= t(matrix(fitmom$G, nrow = J+4)), 
						H 	= fitmom$H,
						Cost = fitmom$f,
						verbose = FALSE,
						timeout = timeout
					),
				silent=FALSE
				)
			if('try-error' %in% class(fit) || fit$IsError){
				lammax_hat = norm(Sn,"2")/(1+sqrt(p/n))
				lambda = NULL
			}
			else{
				wJ = fit$X[1:J]
				ow = cumsum(wJ)
				idw = which(ow>=p/(p+1))
				lammax_hat = tJ[idw[1]]
				idk = sapply(1:p, 
							function(k){ 
								which(ow>=k/(p+1))[1]
							}
						)
				lambda = tJ[idk]
			}
		}
		else{
			wJ = fit$X[1:J]
			ow = cumsum(wJ)
			idw = which(ow>=p/(p+1))
			lammax_hat = tJ[idw[1]]
			idk = sapply(1:p, 
						function(k){
							which(ow>=k/(p+1))[1]
						}
					)
			lambda = tJ[idk]
		}
	}
	else if(method=='mpmo'){
		if(is.null(nmoms)){
			nmoms = 7
		}
		else if(nmoms < 0 ){
			stop("Mumber of moments is a nonegative integer!")
		}
		Sn		= X%*%t(X)/n
		eigs	= eigen(Sn)$values

		fitmom	= .Call("_ARGUEMENTS_MPMOM",
					as.numeric(Sn),
					as.numeric(eigs),
					as.integer(c(n, p, nmoms, J, 0))
					)
		tJ = fitmom$tJ

		fit	<- try(
				linp( E 	= fitmom$E,
					F 	= 1,
					G 	= t(matrix(fitmom$G, nrow = J+nmoms)), 
					H 	= fitmom$H,
					Cost = fitmom$f,
					verbose = FALSE,
					timeout = timeout
				),
			silent=FALSE
			)

		if('try-error' %in% class(fit) ||fit$IsError){
			fitmom	= .Call("_ARGUEMENTS_MPMOM",
						as.numeric(Sn),
						as.numeric(eigs),
						as.integer(c(n, p, 4, J, 0))
						)
			tJ = fitmom$tJ

			fit	<- try(
					linp( E 	= fitmom$E,
						F 	= 1,
						G 	= t(matrix(fitmom$G, nrow = J+4)), 
						H 	= fitmom$H,
						Cost = fitmom$f,
						verbose = FALSE,
						timeout = timeout
					),
				silent=FALSE
				)
			if('try-error' %in% class(fit) || fit$IsError){
				lammax_hat = norm(Sn,"2")/(1+sqrt(p/n))
				lambda = NULL
			}
			else{
				wJ = fit$X[1:J]
				ow = cumsum(wJ)
				idw = which(ow>=p/(p+1))
				lammax_hat = tJ[idw[1]]
				idk = sapply(1:p, 
							function(k){ 
								which(ow>=k/(p+1))[1]
							}
						)
				lambda = tJ[idk]
			}
		}
		else{
			wJ = fit$X[1:J]
			ow = cumsum(wJ)
			idw = which(ow>=p/(p+1))
			lammax_hat = tJ[idw[1]]
			idk = sapply(1:p, 
						function(k){ 
							which(ow>=k/(p+1))[1]
						}
					)
			lambda = tJ[idk]
		}
	}
	else{
		lammax_hat = norm(Sn,"2")/(1+sqrt(p/n))
		lambda = NULL
	}

	return(list(lammax = lammax_hat, lamest = lambda))
}