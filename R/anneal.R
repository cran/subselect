anneal<-function(mat, kmin, kmax=kmin, nsol=1, niter=1000,
exclude=NULL, include=NULL, improvement=TRUE, 
setseed = FALSE, cooling=0.05, temp=1,
coolfreq=1, criterion="RM", pcindices=1:kmax, initialsol=c(0)){


# validation of input

        if (!is.matrix(mat)) 
          stop("Data is missing or is not given in matrix form")
        if (dim(mat)[1] != dim(mat)[2]) {
          mat <- var(mat)
          warning("Data must be given as a covariance or correlation matrix.\n It has been assumed that you wanted the covariance matrix of the \n data matrix supplied")}
        if (kmax < kmin) {
          aux<-kmin
          kmin<-kmax
          kmax<-aux
          warning("Argument kmin should precede argument kmax.\n Since the value of kmin exceeded that of kmax, they have been swapped \n")}
        if (!(as.integer(coolfreq) == coolfreq) | (coolfreq < 1)) stop("\n The cooling frequency must be a non-negative integer")
        p<-dim(mat)[[2]]
        if (kmin >= p) {
             kmin<-p-1
             warning("\n The value of kmin requested is equal to or exceeds the number \n of variables. It has been set at p-1")
            }  
        if (kmax >= p) {
             kmax<-p-1
             warning("\n The value of kmax requested is equal to or exceeds the number \n of variables. It has been set at p-1")
            }
        if (sum(duplicated(c(exclude,include))) > 0)
          stop("\n You have requested that the same variable be both included in, and excluded \n from, the subset.")
         labelsrm<-c("RM","Rm","rm","1",1)
         labelsrv<-c("RV","Rv","rv","2",2)
         labelsgcd<-c("GCD","Gcd","gcd","3",3)
         if (sum(criterion == c(labelsrm,labelsrv,labelsgcd)) == 0) stop("criterion requested is not catered for, or has been misspecified")
         

# initializations for the Fortran subroutine
       
        if (setseed == TRUE) {set.seed(2,kind="default")} 
        if (sum(criterion == labelsrm) > 0) criterio<-1
        if (sum(criterion == labelsrv) > 0) criterio<-2
        if (sum(criterion == labelsgcd) > 0) criterio<-3

        nexclude<-length(exclude)
        if (kmax >= p-nexclude) stop("\n Cardinalities requested are too large for the requested number of excluded variables")
        if (nexclude !=0) exclude<-sort(exclude)  
        exc<-c(0,exclude)
        ninclude<-length(include)
        if (kmin <= ninclude) stop("\n Cardinalities requested are too small for the requested number of included variables")
        if (ninclude !=0) include<-sort(include)
        inc<-c(0,include)
        if (cooling<=0 || cooling >=1) stop("\n values of cooling must be between 0 and 1")
        if (length(pcindices) != kmax) esp<-TRUE else {if (pcindices
        != 1:kmax) esp<-TRUE  else esp<-FALSE}

# initializations when initial solutions are specified; checking for
# nature of initialsol (and how to interpret it) and for
# conflicts with the exclude and include parameters (which inital
# solutions must respect)

        if (initialsol == c(0)) {silog <- FALSE}
        else {silog <- TRUE  # initial solutions have been specified by user          

# checking for the presence of variables that are to be forcefully excluded
            if ((nexclude != 0) & (sum(exclude == rep(initialsol,rep(length(exclude),length(as.vector(initialsol))))) !=0)) stop("\n the specified initial solutions contain variables that are to be excluded")

# how to deal with various formats of input (of initial solutions) for
# a single cardinality

          dimsol<-dim(initialsol)
          if (length(dimsol) > 3) {
            stop("\n Can't handle arrays of more than 3 dimensions")}
          if (kmin == kmax) {

# inital solution is a vector: must be repeated if nsol > 1

           if (is.vector(initialsol)) {
            if (length(initialsol) != kmax) 
             {stop("\n The specified initial and final solutions have different cardinalities")}
            else
             {if (nsol > 1) initialsol<-matrix(nrow=nsol,ncol=kmax,rep(initialsol,rep(nsol,kmax)))}
           }

# initial solution is array?

           if (is.array(initialsol)) {

# initialsol is 3-d array
           if (length(dimsol) == 3) {
            if (dim(initialsol)[[3]] > 1) stop("\n The input array of initial solutions must have dimensions as \n (1 or nsol) x k x 1 when a single cardinality is requested")
            else initialsol<-matrix(nrow=dim(initialsol)[[1]],ncol=dim(initialsol)[[2]],initialsol)
            }
           else

# initialsol is matrix: must be of form nsol x k
 
            {if (dim(initialsol)[[2]] != kmax) stop("\n Input matrix of initial solutions must have as many columns as variables in the requested subset")
             else 
              {if  (dim(initialsol)[[1]] != nsol) {
                if (dim(initialsol)[[1]] == 1) {initialsol<-matrix(nrow=nsol,ncol=dim(initialsol)[[2]],rep(initialsol,rep(nsol,kmax))) }
                else {stop("\n The number of initial solutions can only be 1 or nsol (number of final solutions requested)")}}}}
              }}

# how to deal with various formats of input (initialsol) if more than
# one cardinality is requested

             else  # (if kmax > kmin), i.e., more than one cardinality requested

# initial solution cannot be a single vector

              {if (is.vector(initialsol)) stop("\n There must be initial solutions for all cardinalities requested")

               if (is.array(initialsol)) {

# initial solution 3-d array?

                  if (length(dim(initialsol)) == 3) {
                   if ((dim(initialsol)[[3]] != length(kmin:kmax)) | (dim(initialsol)[[2]] != kmax) | ((dim(initialsol)[[1]] != nsol) & (dim(initialsol)[[1]] > 1))) stop("\n The input array of initial solutions must have dimensions as \n (1 or nsol) x kmax x no. of different cardinalities requested")
                   else 
                    { if (dim(initialsol)[[1]] != nsol) initialsol<-array(dim=c(nsol,kmax,length(kmin:kmax)),rep(initialsol,rep(nsol,kmax*length(kmin:kmax))))
                  }}
                else 

# if initalsol is a matrix

                    {
                    if ((dim(initialsol)[[2]] != kmax) | (dim(initialsol)[[1]] != length(kmin:kmax))) stop("\n A matrix of initial solutions for more than one cardinality must be of type no. of cardinalities x kmax")
                    if (nsol > 1) initialsol<-array(dim=c(nsol,kmax,length(kmin:kmax)),rep(t(initialsol),rep(nsol,kmax*length(kmin:kmax))))
}}}
                 if ((ninclude != 0) & (sum(include == rep(initialsol,rep(length(include),length(as.vector(initialsol))))) != nsol*length(kmin:kmax)*length(include))) stop("\n Not all the specified initial solutions contain the variables that are to be included") 
}

#
# initializations for Fortran output
#

        valores<-rep(0.0,length(kmin:kmax)*nsol)    
        vars<-rep(0,nsol*length(kmin:kmax)*kmax)
        bestval<-rep(0.0,length(kmin:kmax))
        bestvar<-rep(0,kmax*length(kmin:kmax))

#
# call to the Fortran subroutine
#

        Fortout<-.Fortran("anneal",as.integer(criterio),as.integer(p),
          as.double(as.vector(mat)),
          as.integer(kmin),as.integer(kmax),as.double(valores),
          as.integer(vars),as.double(bestval),as.integer(bestvar),
          as.integer(nexclude),as.integer(exc),as.integer(ninclude),
          as.integer(inc),as.integer(nsol),as.integer(niter),
          as.logical(improvement),
          as.double(cooling),as.double(temp),
          as.integer(coolfreq),as.integer(length(pcindices)),
          as.integer(pcindices),as.logical(esp),as.logical(silog),
          as.integer(as.vector(initialsol)),PACKAGE="subselect")

# preparing the output

        valores<-matrix(nrow=nsol,ncol=length(kmin:kmax),Fortout[[6]])
        dimnames(valores)<-list(paste("Solution",1:nsol,sep=" "),paste("card.",kmin:kmax,sep=""))
        variaveis<-array(Fortout[[7]],c(nsol,kmax,length(kmin:kmax)))
        dimnames(variaveis)<-list(paste("Solution",1:nsol,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
        bestval<-Fortout[[8]]
        names(bestval)<-paste("Card",kmin:kmax,sep=".")
        bestvar<-t(matrix(nrow=kmax,ncol=length(kmin:kmax),Fortout[[9]]))
        dimnames(bestvar)<-list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))
        output<-list(variaveis,valores,bestval,bestvar)
        names(output)<-c("subsets","values","bestvalues","bestsets")
        output}


