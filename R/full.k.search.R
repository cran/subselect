full.k.search<-function(mat, k, print.all = FALSE, file="")
{
#  error checking

        if (!is.matrix(mat)) stop("Data is missing or is not given in matrix form")
        if (dim(mat)[1] != dim(mat)[2])  {
             mat<-cov(mat) 
             warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied")
        }

# body of function

	p <- dim(mat)[[2]]
	indices <- c((p - k + 1):p, -1)
	j <- k
	rmmax <- rm.coef(mat, indices[1:k])
	gcdmax <- gcd.coef(mat, indices[1:k], pcindices = seq(1:k))
	rvmax <- rv.coef(mat, indices[1:k])
	rmind <- indices[1:k]
	gcdind <- indices[1:k]
	rvind <- indices[1:k]
        if (print.all) cat(file=file, "   rm   ", "\t",  "  gcd   ", "\t", "   rv   ", "\t", "indices","\n")
	while(j < (k + 1)) {
		j <- 1
		rm <- rm.coef(mat, indices[1:k])
		gcd <- gcd.coef(mat, indices[1:k], pcindices = seq(1:k))
		rv <- rv.coef(mat, indices[1:k])
		if(print.all) {
			cat(file=file, rm, "\t",  gcd, "\t", rv, "\t", indices[1:k],"\n")
		}
		if(rm >= rmmax) {
			rmmax <- rm
			rmind <- indices[1:k]
		}
		if(gcd >= gcdmax) {
			gcdmax <- gcd
			gcdind <- indices[1:k]
		}
		if(rv >= rvmax) {
			rvmax <- rv
			rvind <- indices[1:k]
		}
		while(indices[j] == j) {
			j <- j + 1
		}
		indices[j] <- indices[j] - 1
		if(j > 1 && j <= k) {
			indices[1:(j - 1)] <- (indices[j] - (j - 1)):(indices[j] - 1)
		}
	}
	results <- list(rmmax, rmind, gcdmax, gcdind, rvmax, rvind)
	names(results) <- c("rmmax", "rmindices", "gcdmax", "gcdindices", "rvmax", "rvindices")
	results
}
