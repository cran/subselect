genetic<-function(mat, kmin, kmax=kmin, popsize=100, nger=20,
mutate=FALSE, maxclone=5, exclude=NULL, include=NULL,
improvement=TRUE, printfile=FALSE,
iseed=c(sample(0:4095,3),sample(0:2047,1)*2+1), criterion="RM", pcindices=1:kmax){

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

        if (sum(criterion == labelsrm) > 0) criterio<-1
        if (sum(criterion == labelsrv) > 0) criterio<-2
        if (sum(criterion == labelsgcd) > 0) criterio<-3

        nexclude<-length(exclude)
        if (nexclude !=0) exclude<-sort(exclude)  
        exc<-c(0,exclude)
        ninclude<-length(include)
        if (ninclude !=0) include<-sort(include)
        inc<-c(0,include)
        if (length(pcindices) != kmax) esp<-TRUE else {if (pcindices != 1:kmax) esp<-TRUE  else esp<-FALSE}
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
          as.logical(improvement),as.logical(printfile),
          as.integer(iseed),
          as.integer(length(pcindices)),as.integer(pcindices),as.logical(esp),
          as.integer(kabort),
          PACKAGE="subselect")

# preparing the output

        kabort<-Fortout[[24]]
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


