rm.coef<-function(mat, indices)
{

#   Computes the matrix correlation between data matrices and their 
#   regression on a subset of their variables. Expected input is a
#   variance-covariance (or correlation) matrix. 

#  error checking

        if (!is.matrix(mat)) stop("Data is missing or is not given in matrix form")
        if (dim(mat)[1] != dim(mat)[2])  {
             mat<-cov(mat) 
             warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied")
        }
# function to compute the trace of a saure matrix.

   tr<-function(mat){sum(diag(mat))}

# body of function

	sqrt(tr((mat %*% mat)[indices, indices] %*% solve(mat[indices, indices]))/tr(mat))
}
