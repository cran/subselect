rv.coef<-function(mat, indices)
{
#    Computes Escoufier's RV-coefficient for the configuration of
#    points defined by n observations of a set of p variables, and by
#    the regression of all variables on a subset of k variables
#    (given by \code{indices}). 

#  error checking

        if (!is.matrix(mat)) stop("Data is missing or is not given in matrix form")
        if (dim(mat)[1] != dim(mat)[2])  {
             mat<-cov(mat) 
             warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied")
        }
# function to compute the trace of a square matrix

       tr<-function(mat){sum(diag(mat))}   

# body of function

	mat2 <- (mat %*% mat)[indices, indices]
	invmatk <- solve(mat[indices, indices])
	sqrt(tr(mat2 %*% invmatk %*% mat2 %*% invmatk)/tr(mat %*% mat))
}
