validmat<-function(mat,p,tolval,tolsym)
{
#   Function validmat
#   Checks the main (covariance/total sums of squares and products)
#   matrix for size, symmetry, positive definiteness, etc.
#   To be called from other functions, not directly by the user.

# valid input?

  if (!is.matrix(mat)) {
         stop("Data is missing or is not given in matrix form")}
     if (dim(mat)[1] != dim(mat)[2]) {
         mat<-cor(mat)
         warning("Data must be given as a covariance or correlation matrix. \n It has been assumed that you wanted the correlation matrix of the \n data matrix which was supplied.")
       }

# checking for symmetry in the 'total' input matrix

         maxabssym <- max(abs(mat-t(mat)))
         if (maxabssym > tolsym)
           {stop("\n The covariance/total matrix supplied is not symmetric.\n  Symmetric entries differ by up to ",maxabssym,".")}
         else if (maxabssym > .Machine$double.eps) {
              mat <-(mat+t(mat))/2
                    warning("\n The covariance/total matrix supplied was slightly asymmetric: \n symmetric entries differed by up to ",maxabssym,".\n (less than the 'tolsym' parameter).\n It has been replaced by its symmetric part.\n")
}

  # Positive definiteness, etc.
  
        if (qr(mat)$rank != p) 
            stop("\n The covariance/correlation matrix supplied is not of full rank") 
 	if ( eigen(mat,symmetric=TRUE,only.values=TRUE)$values[p] < -1*tolval )
            stop("\n The covariance/correlation matrix supplied is not positive definite")

# checking for ill-conditioned 'total' matrix

         eigvals<-eigen(mat,symmetric=TRUE)
         if (eigvals$values[p]/eigvals$values[1] < tolval) stop(paste("\n The covariance/correlation matrix supplied has reciprocal condition number \n smaller than the specified threshold of",tolval,". \n Setting a lower value of the 'tolval' function argument may force a solution, \n but numerical accuracy  may be compromised. \n Try a new function call after excluding  variables responsible for the \n (real or approximate) linear dependences. See help(trim.matrix) for  \n assistance in this respect."))
  

}

