trim.matrix <- function(mat,tolval=10*.Machine$double.eps)
{
  p <- dim(mat)[2]
  matindices <- 1:p                   # Pedro 1-07-2013
  mat.eig <- eigen(mat,symmetric=TRUE)
  discard <- rep(FALSE,p)
  newmat <- mat
  newmatindices <- matindices         # Pedro 1-07-2013  
  while(mat.eig$values[p]/mat.eig$values[1] < tolval)
  {
    int <- as.numeric(newmatindices[order(abs(mat.eig$vectors[,p]),decreasing=TRUE)[1]])   
    discard[int] <- TRUE
    newmat <- mat[!discard,!discard]
    newmatindices <- matindices[!discard]    # Pedro 1-07-2013
    p <- p-1
    mat.eig <- eigen(newmat,symmetric=TRUE)   # Pedro 1-07-2013
  }
  size <- dim(newmat)[2]
#  colnames(newmat) <- stored.names[!discard]
  output <- list(newmat,as.numeric(matindices[discard]),colnames(mat)[discard],size)   # Pedro 1-07-2013   
  names(output) <- c("trimmedmat","numbers.discarded","names.discarded","size")
  output
}