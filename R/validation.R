validation<-function(mat, kmin, kmax, exclude, include, criterion, pcindices){

##########################################################
#  general validation of input for all search functions  #
##########################################################


####################################################################
# checking for an input matrix that must be square, of full rank,  #
#    symmetric, positive definite                                  #
####################################################################

         if (!is.matrix(mat)) 
          stop("Data is missing or is not given in matrix form\n") 
         p<-dim(mat)[2] 
         if (qr(mat)$rank != p) stop("\n The covariance/correlation matrix supplied is not of full rank\n")
         if (dim(mat)[1] != p) {
         mat <- var(mat)
         warning("Data must be given as a covariance or correlation matrix.\n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied\n")}
         if ( max( abs(mat-t(mat)) ) > 1E-6 ) 
	    stop("\n The covariance/correlation matrix supplied is not symmetric\n")
  	 if ( eigen(mat,only.values=TRUE)$values[p] < 0 )
            stop("\n The covariance/correlation matrix supplied is not positive definite\n")

#################################################
# checking acceptability of criterion requested #
#################################################
         
         labelsrm<-c("RM","Rm","rm","1",1)
         labelsrv<-c("RV","Rv","rv","2",2)
         labelsgcd<-c("GCD","Gcd","gcd","3",3)
         if (sum(criterion == c(labelsrm,labelsrv,labelsgcd)) == 0) stop("criterion requested is not catered for, or has been misspecified\n")

         if (sum(criterion == labelsrm) > 0) criterio<-1
         if (sum(criterion == labelsrv) > 0) criterio<-2
         if (sum(criterion == labelsgcd) > 0) criterio<-3


         
######################################################################
#  checking consistency of the requested values of kmin, kmax        #
#  (extreme sizes of variable subsets). The validation that kmax<=p  #
#   is made later, considering the no. of excluded variables         #
######################################################################         

         if (!is.numeric(kmin) || !is.numeric(kmax)) stop("\n Arguments kmin and kmax must be numeric.\n Perhaps unnamed arguments in the wrong order?")
         if (kmax < kmin) {
          aux<-kmin
          kmin<-kmax
          kmax<-aux
          warning("the argument kmin should precede the argument kmax.\n Since the value of kmin exceeded that of kmax, they have been swapped \n")}
         if (kmin >= p) {
             kmin<-p-1
             warning("\n The value of kmin requested is equal to or exceeds the number \n of variables. It has been set at p-1. \n")
            }
         if (kmax >= p) {
                         kmax<-p-1
                         warning("\n The value of kmax requested is equal to or exceeds the number \n of variables. It has been set at p-1. \n")
            }

#############################################################
# checking for consistency of requests to exclude/include   #
#  certain variables.                                       #
#############################################################

         if (sum(duplicated(c(exclude,include))) > 0)
          stop("\n You have requested that the same variable be both included in, and excluded \n from, the subset.\n")
         nexclude<-length(exclude)
         if (nexclude !=0) {
         if (kmax >= p-nexclude) {
                                  kmax<-p-nexclude-1
                                  warning("\n Cardinalities requested are too large for the requested number of excluded \n variables, and kmax has been set at p-nexclude-1 \n")}
         exclude<-sort(exclude)}  
         exc<-c(0,exclude)
         ninclude<-length(include)
         if (ninclude !=0){
         if (kmin <= ninclude) {kmin<-ninclude+1
                                warning("\n Cardinalities requested are too small for the requested number of included \n variables, and kmin has been set to ninclude + 1 \n")}
         include<-sort(include)}
         inc<-c(0,include)
         if (kmax<kmin) stop("\n After trying to adapt to the requests for exclusion and inclusion of variables, \n kmax is now smaller than kmin. \n There must be a mistake\n")

######################
# Checking pcindices #
######################


         
# Value always assigned to esp and to non-numeric pcindices to enable passing to parent.frame and to Fortran or C++ routines, for criteria other than GCD

         esp<-FALSE
         if (criterio == 3) { 
                            if (is.null(pcindices)) stop("\n For criterion GCD, argument pcindices must be explicitely set in the leaps \n function, must be non-NULL in other search functions. \n")
                            if (is.numeric(pcindices))  {esp<-TRUE
                                    if  (sum(!(as.integer(pcindices) == pcindices)) > 0) stop("\n The PC indices must be integers.\n")
                                    if (max(pcindices)  >  p) stop("\n PCs of rank larger than the data set were requested. \n")
                                                        }
                            else {if (pcindices != "first_k")
                                    {stop("\n unrecognized value for 'pcindices' argument \n")}}}
         if ((pcindices == "first_k") || is.null(pcindices)) {pcindices<-1:kmax}

#################################
# assigning any changed values  #
#################################
         
         assign("mat",mat,pos=parent.frame())         
         assign("kmax",kmax,pos=parent.frame())         
         assign("kmin",kmin,pos=parent.frame())         
         assign("exclude",exclude,pos=parent.frame())         
         assign("include",include,pos=parent.frame())         
         assign("nexclude",nexclude,pos=parent.frame())         
         assign("ninclude",ninclude,pos=parent.frame())         
         assign("criterio",criterio,pos=parent.frame())    
         assign("pcindices",pcindices,pos=parent.frame())         
         assign("p",p,pos=parent.frame())                  
         assign("exc",exc,pos=parent.frame())         
         assign("inc",inc,pos=parent.frame())         
         assign("esp",esp,pos=parent.frame())
       }


