genetic<-function(mat, kmin, kmax=kmin, popsize=100, nger=100,
mutate=FALSE, mutprob=0.01, maxclone=5, exclude=NULL, include=NULL,
improvement=TRUE, setseed= FALSE,  criterion="RM", pcindices=1:kmax,
initialpop=c(0)){


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
          warning("the argument kmin should precede the argument kmax.\n Since the value of kmin exceeded that of kmax, they have been swapped \n")}
        p<-dim(mat)[[2]]
        if (qr(mat)$rank != p) stop("\n The covariance/correlation matrix supplied is not of full rank") 
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
        if ((mutprob < 0) | (mutprob > 1)) stop("\n The mutation probability parameter (mutprob) must be between 0 and 1")

        if (setseed == TRUE) set.seed(2,kind="default")
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
        if (length(pcindices) != kmax) esp<-TRUE else {if (pcindices != 1:kmax) esp<-TRUE  else esp<-FALSE}


# initializations when initial population has been specified; checking for
# nature of initialpop (and how to interpret it) and for
# conflicts with the exclude and include parameters (which initial
# population must respect)

        if (initialpop == c(0)) {pilog <- FALSE}
        else {pilog <- TRUE  # initial population has been specified by user          

# checking for the presence of variables that are to be forcefully
# excluded

            if ((nexclude != 0) & (sum(exclude == rep(initialpop,rep(length(exclude),length(as.vector(initialpop))))) !=0)) stop("\n the specified initial population contains variables that are to be excluded")

# how to deal with various formats of input (of initial population) for
# a single cardinality

            dimpop<-dim(initialpop)
            if (length(dimpop) > 3) stop("\n Can't handle arrays of more than 3 dimensions")
            if (kmin == kmax) {

# inital solution is a vector: not acceptable

               if (is.vector(initialpop)) {
                     {stop("\n The specified initial population must have different k-subsets")}
                    }

# initial solution is array?

              if (is.array(initialpop)) {

# initialpop is 3-d array
              if (length(dimpop) == 3) {
                 if (dimpop[[3]] > 1) stop("\n The input array of initial population must have dimensions as \n popsize x k x 1 when a single cardinality is requested")
                 else initialpop<-matrix(nrow=dimpop[[1]],ncol=dimpop[[2]],initialpop)
              }
              else # initialpop is matrix: must be of form popsize x k
 
                 {if (dimpop[[2]] != kmax) stop("\n Input matrix of initial solutions must have as many columns as variables in the requested subset")
                  else 
                      {if  (dimpop[[1]] != popsize) {
                        stop("\n The number of initial solutions can only be 1 or popsize (number of final solutions requested)")}}}
              }}

# how to deal with various formats of input (initialpop) if more than
# one cardinality is requested

             else  #(if kmax > kmin), i.e., more than one cardinality requested

# initial population must be 3-d array

              {if (length(dimpop) < 3) stop("\n There must be initial populations for all cardinalities requested")
               else 
               {if ((dimpop[[3]] != length(kmin:kmax)) | (dimpop[[2]] != kmax) | (dimpop[[1]] != popsize)) stop("\n The input array of initial solutions must have dimensions as \n popsize x kmax x no. of different cardinalities requested")
                  }
                 if ((ninclude != 0) & (sum(include == rep(initialpop,rep(length(include),length(as.vector(initialpop))))) != popsize*length(kmin:kmax)*length(include))) stop("\n Not all the specified initial solutions contain the variables that are to be included") 
}}


# initializations for the Fortran subroutine


        valores<-rep(0.0,length(kmin:kmax)*popsize)
        vars<-rep(0,popsize*length(kmin:kmax)*kmax)
        bestval<-rep(0.0,length(kmin:kmax))
        bestvar<-rep(0,kmax*length(kmin:kmax))
        kabort<-kmax+1

# call to the Fortran subroutine

        Fortout<-.Fortran("genetic",as.integer(criterio),as.integer(p),
          as.double(as.vector(mat)),
          as.integer(kmin),as.integer(kmax),as.double(valores),
          as.integer(vars),as.double(bestval),as.integer(bestvar),
          as.integer(nexclude),as.integer(exc),as.integer(ninclude),
          as.integer(inc),as.integer(popsize),as.integer(nger),
          as.integer(maxclone),as.logical(mutate),
          as.double(mutprob),as.logical(improvement),
          as.integer(length(pcindices)),as.integer(pcindices),as.logical(esp),
          as.integer(kabort),as.logical(pilog),
          as.integer(as.vector(initialpop)),
          PACKAGE="subselect")

# preparing the output

        kabort<-Fortout[[23]]
        valores<-matrix(nrow=popsize,ncol=length(kmin:kmax),Fortout[[6]])
        dimnames(valores)<-list(paste("Solution",1:popsize,sep=" "),paste("card.",kmin:kmax,sep=""))
        variaveis<-array(Fortout[[7]],c(popsize,kmax,length(kmin:kmax)))
        dimnames(variaveis)<-list(paste("Solution",1:popsize,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
        bestval<-Fortout[[8]]
        names(bestval)<-paste("Card",kmin:kmax,sep=".")
        bestvar<-t(matrix(nrow=kmax,ncol=length(kmin:kmax),Fortout[[9]]))
        dimnames(bestvar)<-list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))
        output<-list(variaveis[,1:(kabort-1),1:(kabort-kmin)],valores[,1:(kabort-kmin)],bestval[1:(kabort-kmin)],bestvar[1:(kabort-kmin),1:(kabort-1)])
        names(output)<-c("subsets","values","bestvalues","bestsets")
        if (kabort > kmin) output}


