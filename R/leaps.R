leaps <- function(mat,kmin=1,kmax=ncol(mat)-1,nsol=1,exclude=NULL,include=NULL,
                  criterion="RM",pcindices="first_k",timelimit=15, tolval=10*.Machine$double.eps)
{

#############################################################
# validation specific for the input to the leaps function  #
#############################################################

	if (timelimit <= 0) {stop("\n The time limit argument must be a positive real number")} 


###############################
# general validation of input #
###############################

        validation(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval)



##########################################
# Initializations for the C++ subroutine #
##########################################

	Si = solve(mat)
	if (criterio == 2)  S2 = mat %*% mat
	else S2 = NULL
	if (criterio == 3)  {
	  Spectd = eigen(mat,TRUE)
	  Segval = Spectd$values[pcindices]
	  Segvct = Spectd$vectors[,pcindices]
        }
        else  {
           Segval = NULL
           Segvct = NULL
        }

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
            as.double (S2),
            as.double (Si),
            as.double (Segval),
            as.double (Segvct),
            as.integer(kmin),
            as.integer(kmax),
            as.integer(nsol),
            as.integer(exclude),
            as.integer(include),
            as.integer(nexclude),
            as.integer(ninclude),
            as.character(criterion),
            as.logical(esp),
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
         output <- c(Cout[20:23],match.call())
         names(output) <- c("subsets","values","bestvalues","bestsets","call")
         output

}


