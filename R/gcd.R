gcd.coef<-function(mat, indices, pcindices = seq(1:length(indices)))
{
#   
#       calcula o GCD entre um subconjunto ("indices") de variaveis
#       e um subconjunto ("pcindices") das CPs de todas as variaveis, 
#       cuja matriz de covariancias e "mat".

#  error checking

        if (!is.matrix(mat)) stop("Data is missing or is not given in matrix form")
        if (dim(mat)[1] != dim(mat)[2])  {
             mat<-cov(mat) 
             warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied")
        }
# auxiliary function to compute a low-rank approximation to a matrix,
# based on its Singular Value Decomposition

svdapprox<-function(mat, indices)
{
	dvsmat <- svd(mat)
	t(dvsmat$v[, indices] %*% (t(dvsmat$u[, indices]) * dvsmat$d[indices]))
}


  
# body of function

	sum(diag(solve(mat[indices, indices]) %*% svdapprox(mat, pcindices)[indices, indices]))/sqrt(
length(indices) * length(pcindices))
}
