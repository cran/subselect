leaps <- function(mat,kmin=1,kmax=ncol(mat)-1,nsol=1,exclude=NULL,include=NULL,
                  criterion="RM",pcindices=NULL,timelimit=15)
{

#############################################################
# validation specific for the input to the leaps function  #
#############################################################

 if ((criterion == "3") || (criterion == "GCD") || (criterion == "gcd") || (criterion == "Gcd"  || (criterion == 3))) {  
        if (!is.null(pcindices)) {if (sum(pcindices == "first_k") > 0) stop("\n The 'first_k' option is not available in the 'leaps' function: \n PC indices must be explicitely set. \n")}}
	if (timelimit <= 0) {stop("\n The time limit argument must be a positive real number")} 


###############################
# general validation of input #
###############################

        validation(mat, kmin, kmax, exclude, include, criterion, pcindices)



##########################################
# Initializations for the C++ subroutine #
##########################################

         klength    <- kmax-kmin+1

         subsets    <- integer(nsol*kmax*klength);   dim(subsets)    <- c(nsol,kmax,klength)
         values     <- double(nsol*klength);         dim(values)     <- c(nsol,klength)
         bestvalues <- double(klength);             
         bestsets   <- integer(klength*kmax);        dim(bestsets)   <- c(klength,kmax) 

         dimnames(subsets)   <-  list(paste("Solution",1:nsol,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
         dimnames(values)    <-  list(paste("Solution",1:nsol,sep=" "),paste("card.",kmin:kmax,sep=""))
         names(bestvalues)   <-  paste("Card",kmin:kmax,sep=".")
         dimnames(bestsets)  <-  list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))


############################
# Call to the C subroutine #
############################

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

#######################################
# Preparing and returning the output  #
#######################################
        
	 if (Cout$found == FALSE) {
	    stop("\n Leaps was not able to complete the search within the specified time limit.\n Either increase this limit or try one of the available meta-heuristics \n") }
         output <- c(Cout[15:18],match.call())
         names(output) <- c("subsets","values","bestvalues","bestsets","call")
         output

}


