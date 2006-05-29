genetic<-function(mat, kmin, kmax=kmin, popsize=100, nger=100,
mutate=FALSE, mutprob=0.01, maxclone=5, exclude=NULL, include=NULL,
improvement=TRUE, setseed= FALSE,  criterion="RM", pcindices="first_k",
initialpop=NULL, force=FALSE, H=NULL, r=0,tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps){


        
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
        maxnovar = 400
        if ((p > maxnovar) & (force==FALSE)) stop("\n For very large data sets, memory problems may crash the R session. \n To proceed anyways, repeat the function call with \n the argument 'force' set to 'TRUE' (after saving anything important \n from the current session)\n")


######################################################################
# Parameter validation if the criterion is one of "TAU_2", "XI_2",   #
# "ZETA_2" or "CCR1_2"                                               #
######################################################################

if (criterion == "TAU_2" || criterion == "XI_2" || criterion ==
"ZETA_2" || criterion == "CCR1_2") validnovcrit(mat,criterion,H,r,p,tolval,tolsym)


##############################################################
# more specific validation of input for the genetic function #
##############################################################

        validgenetic(kmin, kmax, popsize, mutprob, exclude, nexclude, include, ninclude, initialpop)


##############################################
# initializations for the Fortran subroutine #
##############################################

        if (setseed == TRUE) set.seed(2,kind="default")
        valores<-rep(0.0,length(kmin:kmax)*popsize)
        vars<-rep(0,popsize*length(kmin:kmax)*kmax)
        bestval<-rep(0.0,length(kmin:kmax))
        bestvar<-rep(0,kmax*length(kmin:kmax))
        kabort<-kmax+1
        if (criterio == 3) {
          decespectral<-eigen(mat,symmetric=TRUE)
          valp<-decespectral$values
          vecp<-decespectral$vectors
               }
           else {
              valp<-rep(0,p)
              vecp<-matrix(nrow=p,ncol=p,rep(0,p*p))
                 }


##################################
# call to the Fortran subroutine #
##################################

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
          as.integer(as.vector(initialpop)),as.double(valp),
          as.double(as.vector(vecp)),as.double(as.vector(H)),as.integer(r),PACKAGE="subselect")

########################
# preparing the output #
########################

        kabort<-Fortout[[23]]
        valores<-matrix(nrow=popsize,ncol=length(kmin:kmax),Fortout[[6]])
        dimnames(valores)<-list(paste("Solution",1:popsize,sep=" "),paste("card.",kmin:kmax,sep=""))
        variaveis<-array(Fortout[[7]],c(popsize,kmax,length(kmin:kmax)))
        dimnames(variaveis)<-list(paste("Solution",1:popsize,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
        bestval<-Fortout[[8]]
        names(bestval)<-paste("Card",kmin:kmax,sep=".")
        bestvar<-t(matrix(nrow=kmax,ncol=length(kmin:kmax),Fortout[[9]]))
        dimnames(bestvar)<-list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))
        output<-list(variaveis[,1:(kabort-1),1:(kabort-kmin)],valores[,1:(kabort-kmin)],bestval[1:(kabort-kmin)],bestvar[1:(kabort-kmin),1:(kabort-1)],match.call())
        names(output)<-c("subsets","values","bestvalues","bestsets","call")
        if (kabort > kmin) output}

