trim.matrix<-function(mat, tolval=10*.Machine$double.eps){

        if (tolval < 0) stop("\n The 'tolval' argument must be non-negative.")

        p<-dim(mat)[2]
        if ((dim(mat)[1] != p) || (sum(t(mat) != mat) >0)) stop("\n Data must be given as a covariance or correlation matrix.\n The input matrix is rectangular or not symmetric.\n")
        stored.names<-colnames(mat)
        colnames(mat)<-1:p
        mat.eig<-eigen(mat,symmetric=TRUE)
        discard<-rep(FALSE,p)
        newmat<-mat
        while (mat.eig$values[p]/mat.eig$values[1] < tolval) {
              int<-as.numeric(colnames(newmat)[order(abs(mat.eig$vectors[,p]),decreasing=TRUE)[1]])
              discard[int]<-TRUE
              newmat<-mat[!discard,!discard]
              p<-p-1
              mat.eig<-eigen(newmat)
                     }
       size<-dim(newmat)[2]       
       colnames(newmat)<-stored.names[!discard] 
       output<-list(newmat,as.numeric(colnames(mat)[discard]),stored.names[discard],size)
       names(output)<-c("trimmedmat","numbers.discarded","names.discarded","size")
       output}
