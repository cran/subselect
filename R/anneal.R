anneal<-function(mat, kmin, kmax=kmin, nsol=1, niter=1000,
exclude=NULL, include=NULL, improvement=TRUE, 
setseed = FALSE, cooling=0.05, temp=1,
coolfreq=1, criterion="RM", pcindices="first_k", initialsol=NULL){

###############################
# general validation of input #
###############################

        validation(mat, kmin, kmax, exclude, include, criterion, pcindices)

##########################################################################
# more specific validation of input for the anneal and improve functions #
##########################################################################

        implog<-FALSE
        validannimp(kmin, kmax, nsol, exclude, nexclude, include, ninclude, initialsol,implog)


#############################################################
# validation specific for the input to the anneal function  #
#############################################################

        if (!(as.integer(coolfreq) == coolfreq) | (coolfreq < 1)) stop("\n The cooling frequency must be a non-negative integer")
        if (cooling<=0 || cooling >=1) stop("\n values of cooling must be between 0 and 1")



##########################################
# initializations for Fortran subroutine #
##########################################

        if (setseed == TRUE) {set.seed(2,kind="default")} 
        valores<-rep(0.0,length(kmin:kmax)*nsol)    
        vars<-rep(0,nsol*length(kmin:kmax)*kmax)
        bestval<-rep(0.0,length(kmin:kmax))
        bestvar<-rep(0,kmax*length(kmin:kmax))

##################################
# call to the Fortran subroutine #
##################################

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

########################
# preparing the output #
########################

        valores<-matrix(nrow=nsol,ncol=length(kmin:kmax),Fortout[[6]])
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


