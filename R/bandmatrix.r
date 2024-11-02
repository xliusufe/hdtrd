bandmatrix <- function(rho, p, T = 5){

	if(is.null(rho)){
		stop("rho is NULL !")
	}
	if(is.null(p)){
		stop("p is NULL !")
	}

	sighalf = .Call("BANDMATRIX_",
				as.numeric(rho),
				as.integer(c(p, T))
			)
	sighalf = matrix(sighalf, ncol = p)
	sigma 	= t(sighalf)%*%sighalf

	return(list(sighalf	= sighalf,
				sigma	= sigma
				)
			)
}