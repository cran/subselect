genetic<-function(mat, kmin, kmax=kmin, popsize=100, nger=100,
mutate=FALSE, mutprob=0.01, maxclone=5, exclude=NULL, include=NULL,
improvement=TRUE, setseed= FALSE,  criterion="default", pcindices="first_k",
initialpop=NULL, force=FALSE, H=NULL, r=0,tolval=1000*.Machine$double.eps,tolsym=1000*.Machine$double.eps){


#########################################################################################
# auxiliary  variables (includes declarations to avoid "no visible binding" NOTE)       #
#########################################################################################

	p <- ncol(mat)    				# Number of original variables
	nexclude <- length(exclude)     		# Number of excluded variables
	ninclude <- length(include)     		# Number of included variables				 
	if (pcindices!="first_k") esp <- TRUE		# The user has specified the set of Principal Components to be used with the GCD criterion
	else esp <- FALSE				# The user has not specified the set of Principal Components to be used with the GCD criterion
	if (!is.null(initialpop)) pilog <- TRUE		# The user has specified an initial population
	else pilog <- FALSE				# The user has not specified an initial population
	exc <- exclude					# Name for vector of excluded variables after 'validation' routine 
	inc <- include					# Name for vector of included variables after 'validation' routine 

        
###############################
# general validation of input #
###############################

  
        initialization(mat, criterion, r)
        validation(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval,tolsym)
        maxnovar = 400
        if ((p > maxnovar) & (force==FALSE)) stop("\n For very large data sets, memory problems may crash the R session. \n To proceed anyways, repeat the function call with \n the argument 'force' set to 'TRUE' (after saving anything important \n from the current session)\n")


######################################################################
# Parameter validation if the criterion is one of "TAU_2", "XI_2",   #
# "ZETA_2" or "CCR1_2"  or "WALD"                                    #
######################################################################

if (criterion == "TAU_2" || criterion == "XI_2" || criterion ==
"ZETA_2" || criterion == "CCR1_2" || criterion == "WALD") validnovcrit(mat,criterion,H,r,p,tolval,tolsym)


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

	if (criterion == "WALD")   {

###       Convert a min Wald problem into an (artificial) equivalent max "XI2" problem
        
            Waldval <- sum(diag(solve(mat,H)))
                H <- H / Waldval 
		criterion <- "XI_2"
		criterio <- 5
		Walddec <- TRUE
	}
        else Walddec <- FALSE

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
        output<-list(variaveis[,1:(kabort-1),1:(kabort-kmin),drop=FALSE],valores[,1:(kabort-kmin),drop=FALSE],bestval[1:(kabort-kmin)],bestvar[1:(kabort-kmin),1:(kabort-1),drop=FALSE],match.call())
        names(output)<-c("subsets","values","bestvalues","bestsets","call")
        if (Walddec)  {
		criterion <- "WALD"
		criterio <- 8
		validvalues <- output$values[output$values > 0.]
		output$values[output$values > 0.] <- rep(Waldval,length(validvalues)) - validvalues * Waldval
		output$bestvalues <- rep(Waldval,kmax-kmin+1) - output$bestvalues * Waldval
        }
    if (kabort > kmin){
        output
       } else
       {
        output  <- NA
        return(invisible(output))
       }
      }


