ccr12.coef<-function(mat, H=NULL,r=0,indices,tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps)
{
#   Function ccr1_2.coef
#   Computes the first squared canonical correlation.  
#   This criterion is equivalent to the maximization of Roy first root.
#   Expected input: a variance-covariance (or correlation) matrix, 
#   Effect descrption matrix (H) and its rank (r)

#  error checking

#  mat and indices

  if  (sum(!(as.integer(indices) == indices)) > 0) stop("\n The variable indices must be integers")

  p <- dim(mat)[2]
  validmat(mat,p,tolval,tolsym)  


# checks on r and H

  validnovcrit(mat,criterion="CCR1_2",H,r,p,tolval,tolsym)
  
#  Computing the criterion value  

      rm.1d<-function(mat,H,r,indices){
		
	Re(eigen(H[indices,indices]%*%solve(mat[indices,indices]))$values[1])
		
       }
      dimension<-length(dim(indices))
      if (dimension > 1){
         rm.2d<-function(mat,H,r,subsets){
             apply(subsets,1,function(indices){
			
			rm.1d(mat,H,r,indices)})
            }  
             if (dimension > 2) {
               rm.3d<-function(mat,H,r,array3d){
                   apply(array3d,3,function(subsets){rm.2d(mat,H,r,subsets)})
                 }
               output<-rm.3d(mat,H=H,r,indices)
              }
              if (dimension == 2) {output<-rm.2d(mat=mat,H,r,indices)}
      }

      if (dimension < 2) {output<-rm.1d(mat,H,r,indices)}
      output
}
