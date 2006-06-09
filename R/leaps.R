leaps<-function(mat,kmin=1,kmax=ncol(mat)-1,nsol=1,exclude=NULL,include=NULL,criterion="default",pcindices="first_k",timelimit=15,H=NULL,r=0,tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps)
{


        
#####################################
#  set parameters to default values #
#####################################

	if (r==0 && criterion=="default")  criterion <- "RM"
        if (r>0  && criterion=="default")  criterion <- "TAU_2"
        p <- nrow(mat)

        
###############################
# general validation of input #
###############################


        validation(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval,tolsym)



#############################################################
# validation specific for the input to the leaps function  #
#############################################################

	if (timelimit <= 0) stop("\n The time limit argument must be a positive real number\n")
        if ((criterion =="CCR1_2") && (r > 3)) stop("\n The 'leaps' function does not accept, for the CCR1_2 criterion, \n an effects matrix (H) with  rank greater than 3\n")

        
######################################################################
# Parameter validation if the criterion is one of "TAU_2", "XI_2",   #
# "ZETA_2" or "CCR1_2"                                               #
######################################################################

if (criterion == "TAU_2" || criterion == "XI_2" || criterion ==
"ZETA_2" || criterion == "CCR1_2") validnovcrit(mat,criterion,H,r,p,tolval,tolsym)


##########################################
# Initializations for the C++ subroutine #
##########################################

       # Normalize matrices in order to improve numerical accuracy

        nfactor <- 1./max(mat)
	mat <- nfactor * mat
        if (r > 0) H <- nfactor * H  
	
	# Matrix initializations
	
	Si <- solve(mat)

	if (criterion == "RV")  S2 <- mat %*% mat
	else S2 <- NULL

	if (criterion == "GCD")  {
	  SSpectd <- eigen(mat,symmetric=TRUE)
	  Segval <- SSpectd$values[pcindices]
	  Segvct <- SSpectd$vectors[,pcindices]
        }
        else  {
           Segval <- NULL
           Segvct <- NULL
        }
 
	if (criterion == "TAU_2" || criterion == "ZETA_2" || (criterion == "CCR1_2" && r > 1) ) {
	   E <- mat - H
           Ei <- solve(E)
        }
	else  E <- Ei <- NULL
 
	if (criterion == "XI_2" || criterion == "ZETA_2" || criterion == "CCR1_2") 	{
	  HSpectd <- eigen(H,symmetric=TRUE)
	  if (r > 1) Hegvct <- HSpectd$vectors[,1:r] %*% sqrt(diag(HSpectd$values[1:r])) 
	  else Hegvct <- HSpectd$vectors[,1] * sqrt(HSpectd$values[1]) 
        }
	else  Hegvct <- NULL
	if (criterion == "XI_2" || criterion == "CCR1_2") HegvctTinv <- Si %*% Hegvct
	else HegvctTinv <- NULL
	if (criterion == "ZETA_2" || (criterion == "CCR1_2" && r == 3)) HegvctEinv <- Ei %*% Hegvct  
	else HegvctEinv <- NULL

	if ( criterion == "TAU_2" || (criterion == "CCR1_2" && r > 1) ) Wilksval <- det(E) / det(mat)
	else Wilksval <- 0.
	if ( criterion == "XI_2" || criterion == "CCR1_2"  ) HSi <- H %*% Si	
	if ( criterion == "XI_2" || (criterion == "CCR1_2" && r > 1) ) BartPival <- sum(diag(HSi))	
	else BartPival <- 0.
	if ( criterion == "ZETA_2" || (criterion == "CCR1_2" && r == 3) ) LawHotval <- sum(diag(H %*% Ei))	
	else LawHotval <- 0.
	if ( criterion == "CCR1_2") CCR12val <- as.numeric(eigen(HSi,symmetric=FALSE,only.values=TRUE)$values[1])	
	else CCR12val <- 0.
 
	# Memory initializations
	
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
            as.double(mat),
            as.double(S2),
            as.double(Si),
            as.double(Segval),
            as.double(Segvct),
            as.double(E),
            as.double(Ei),
            as.double(Hegvct),
            as.double(HegvctTinv),
            as.double(HegvctEinv),
	    as.double(Wilksval),			
	    as.double(BartPival),			
	    as.double(LawHotval),			
	    as.double(CCR12val),			
	    as.integer(r),
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
	    warning("\n Leaps was not able to complete the search within the specified time limit.\n Either increase this limit or try one of the available meta-heuristics\n") 
            NULL
         }
	 else {
           output <- c(Cout[30:33],match.call())
           names(output) <- c("subsets","values","bestvalues","bestsets","call")
           output
         } 
}
