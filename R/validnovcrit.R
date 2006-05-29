validnovcrit<-function(mat,criterion,H,r,p,tolval,tolsym){



################################################################
# validation of input that is specific to the criteria "TAU_2",#
#    "XI_2","ZETA_2" and "CCR1_2"                              #    
################################################################

       if (!is.numeric(r)) stop("\n Argument r must be numeric.\n Perhaps unnamed arguments in the wrong order?\n")	
             
       if (r<=0) stop("\n Argument r must be a positive integer (expected rank of effects matrix H)\n")
       if (r!=floor(r)) {
                   warning("\n The specified expected rank 'r' of the effects matrix must be an integer. \n The specified value of the parameter 'r' has been truncated.\n")
                   r <- floor(r)
                        }
     
      if (!is.matrix(H)) stop("\nEffect description data (H) is missing or is not given in matrix form \n")
      if (nrow(H) != p || ncol(H) != p) stop("\nEffect description matrix (H) does not have correct dimensions \n")
                
# checking symmetry of the H effects matrix
                
                maxabssym <- max(abs(H-t(H)))
                if (maxabssym > tolsym)
                      {stop("\n The effect description matrix (H) supplied is not symmetric.\n Symmetric entries differ by up to ",maxabssym,".")}
                else if (maxabssym > .Machine$double.eps) {
                    H <-(H+t(H))/2
                    warning("\n The effect description matrix (H) supplied was slightly asymmetric: \n symmetric entries differed by up to ",maxabssym,".\n (less than the 'tolsym' parameter).\n The H matrix has been replaced by its symmetric part.\n")
                }

# checking whether H matrix has the specified rank
                
		if (qr(H)$rank != r) {
                  if (criterion == "CCR1_2") stop("\nEffect description matrix (H) does not have the specified rank\n")
                  else warning("\n The expected rank of the effect matrix supplied by the user is ",r," but \n the matrix rank (as computed by R's qr function) is ",qr(H)$rank,".\n")}

# A VER - tirou-se a ref a XI2 para tentar evitar valore absurdos de E. Mas para "E"
# plausivel, poderia ser XI2. Tentar arrumar isto.
#		if (criterion != "XI_2" &&  eigen(mat-H)$value[p] < -1*tolval) stop("\nError #matrix (E=T-H) is not of full rank (a singular E is not acceptable for this criterion)\n")

		if (eigen(mat-H)$value[p] < -1*tolval) stop("\nError matrix (E=T-H) is not of full rank (a singular E is not acceptable for\n this criterion). Check whether the efects (H) matrix has been correctly specified.\n")
	
}

