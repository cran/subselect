improve<-function(mat, kmin, kmax=kmin, nsol=1, exclude=NULL,
include=NULL, setseed = FALSE, criterion="RM", pcindices="first_k",
initialsol=NULL, force=FALSE, tolval=.Machine$double.eps){


###############################
# general validation of input #
###############################

        validation(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval)
        maxnovar = 400
        if ((p > maxnovar) & (force==FALSE)) stop("\n For very large data sets, memory problems may crash the R session. \n To proceed anyways, repeat the function call with \n the argument 'force' set to 'TRUE' (after saving anything important \n from the current session)\n")

##########################################################################
# more specific validation of input for the anneal and improve functions #
##########################################################################

        implog<-TRUE
        validannimp(kmin, kmax, nsol, exclude, nexclude, include, ninclude, initialsol, implog)



###########################################
# initializations for Fortran subroutine  #
###########################################

        if (setseed == TRUE) set.seed(2,kind="default")
        valores<-rep(0.0,(kmax-kmin+1)*nsol)
        vars<-rep(0,nsol*length(kmin:kmax)*kmax)
        bestval<-rep(0.0,length(kmin:kmax))
        bestvar<-rep(0,kmax*length(kmin:kmax))
        if (criterio == 3) {
          decespectral<-eigen(mat,symmetric=TRUE)
          valp<-decespectral$values
          vecp<-decespectral$vectors
               }
           else {
              valp<-rep(0,p)
              vecp<-matrix(nrow=p,ncol=p,rep(0,p*p))
                 }

###############################
# call to Fortran subroutine  #
###############################

        Fortout<-.Fortran("improve",as.integer(criterio),as.integer(p),
          as.double(as.vector(mat)),
          as.integer(kmin),as.integer(kmax),as.double(valores),
          as.integer(vars),as.double(bestval),as.integer(bestvar),
          as.integer(nexclude),as.integer(exc),as.integer(ninclude),
          as.integer(inc),as.integer(nsol),
          as.integer(length(pcindices)),as.integer(pcindices),
          as.logical(esp),as.logical(silog),
          as.integer(as.vector(initialsol)),as.double(valp),
          as.double(as.vector(vecp)),PACKAGE="subselect")

########################
# preparing the output #
########################

        valores<-matrix(ncol=length(kmin:kmax),nrow=nsol,Fortout[[6]])
        dimnames(valores)<-list(paste("Solution",1:nsol,sep=" "),paste("card.",kmin:kmax,sep=""))
        variaveis<-array(Fortout[[7]],c(nsol,kmax,length(kmin:kmax)))
        dimnames(variaveis)<-list(paste("Solution",1:nsol,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
        bestval<-Fortout[[8]]
        names(bestval)<-paste("Card",kmin:kmax,sep=".")
        bestvar<-t(matrix(nrow=kmax,ncol=length(kmin:kmax),Fortout[[9]]))
        dimnames(bestvar)<-list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))
        output<-list(variaveis,valores,bestval,bestvar,match.call())
        names(output)<-c("subsets","values","bestvalues","bestsets","call")
        output}



