leaps <- function(mat,kmin=1,kmax=ncol(mat)-1,nsol=1,exclude=NULL,include=NULL,
                  criterion="RM",pcindices=1:kmin,timelimit=15)
{

# Validation of input


        if (!is.matrix(mat)) 
            stop("Data is missing or is not given in matrix form")

   	p <- ncol(mat)
        if (p != nrow(mat)) {
            mat <- var(mat)
            warning("Data must be given as a covariance or correlation matrix.\n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied")
        } 

        if (qr(mat)$rank != p) 
            stop("\n The covariance/correlation matrix supplied is not of full rank") 
  	if ( max( abs(mat-t(mat)) ) > 1E-6 ) 
	    stop("\n The covariance/correlation matrix supplied is not symmetric")
  	if ( eigen(mat,only.values=TRUE)$values[p] < 0 )
            stop("\n The covariance/correlation matrix supplied is not positive definite")

        if (kmax < kmin) {
            aux  <- kmin
            kmin <- kmax
            kmax <- aux
            warning("Argument kmin should precede argument kmax.\n Since the value of kmin exceeded that of kmax, they have been swapped \n")
        }  
        nexclude <- length(exclude)
        if (kmax >= p - nexclude) {
	    kmax <- p - nexclude - 1	
            if (nexclude == 0) 
                warning("\n The value of kmax requested is equal to or exceeds the number\n of variables. It has been set at p-1")
            else 
	        warning("\n The value of kmax requested is too large for the requested number\n of excluded variables. It has been set at p - # excluded variables - 1\n")
        }	 
        if (kmin >= p) {
	    kin <- p - 1	
            warning("\n The value of kmin requested is equal to or exceeds the number\n of variables. It has been set at p-1")
        }	
        ninclude <- length(include)
        if (kmin <= ninclude)  {
	    kmin <- ninclude + 1	
            warning("\n The value of kmin requested is too small for the requested number\n of included variables. It has been set at # included variables + 1\n")
        }
        if (sum(duplicated(c(exclude,include))) > 0) {
            stop("\n You have requested that the same variable be both included in, and excluded\n from, the subset.")}

        labelsrm<-c("RM","Rm","rm","1",1)
        labelsrv<-c("RV","Rv","rv","2",2)
        labelsgcd<-c("GCD","Gcd","gcd","3",3)
        if (sum(criterion == c(labelsrm,labelsrv,labelsgcd)) == 0) {
            stop("\n Criterion requested is not catered for, or has been misspecified")}
        if  (sum(!(as.integer(pcindices) == pcindices)) > 0) {stop("\n The PC indices must be integers")}

	if (timelimit <= 0) {stop("\n The time limit argument must be a positive real number")} 


#
# Initializations for C output
#

         klength    <- kmax-kmin+1

         subsets    <- integer(nsol*kmax*klength);   dim(subsets)    <- c(nsol,kmax,klength)
         values     <- double(nsol*klength);         dim(values)     <- c(nsol,klength)
         bestvalues <- double(klength);             
         bestsets   <- integer(klength*kmax);        dim(bestsets)   <- c(klength,kmax) 

         dimnames(subsets)   <-  list(paste("Solution",1:nsol,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
         dimnames(values)    <-  list(paste("Solution",1:nsol,sep=" "),paste("card.",kmin:kmax,sep=""))
         names(bestvalues)   <-  paste("Card",kmin:kmax,sep=".")
         dimnames(bestsets)  <-  list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))


#
# Call to the C subroutine
#

 	 Cout <- .C("leaps",
            as.double (mat),
            as.integer(kmin),
            as.integer(kmax),
            as.integer(nsol),
            as.integer(exclude),
            as.integer(include),
            as.integer(nexclude),
            as.integer(ninclude),
            as.character(criterion),
            as.integer(pcindices),
            as.integer(length(pcindices)),
            as.integer(p),
	    as.double(timelimit),
	    found = logical(1),	    
            subsets,
            values,
            bestvalues,
            bestsets,
            PACKAGE="subselect"
        ) 


# Preparing and returning the output

        
	 if (Cout$found == FALSE) {
	    stop("\n Leaps was not able to complete the search within the specified time limit.\n Either increase this limit or try one of the available meta-heuristics") }
         if ((Cout[9] == "GCD") || (Cout[9] == "gcd") || (Cout[9] == "Gcd") || (Cout[9] == "3")) {
            warning("\n If no pcindices were specified, the values of the GCD compare \n each k-variable subset with the first kmin PCs \n")}
         output <- Cout[15:18]
         names(output) <- c("subsets","values","bestvalues","bestsets")
         output

}


