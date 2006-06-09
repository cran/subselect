lmHmat <- function(x,...) {
  if (is.null(class(x))) class(x) <- data.class(x)
  UseMethod("lmHmat",x) 
}   

ldaHmat <- function(x,...) {
  if (is.null(class(x))) class(x) <- data.class(x)
  UseMethod("ldaHmat",x)
}

glhHmat <- function(x,...) { 
  if (is.null(class(x))) class(x) <- data.class(x)
  UseMethod("glhHmat",x)
}
  
lmHmat.default <- function(x,y,...)
{
   if  ( !is.matrix(x) || ( !is.matrix(y) && !is.vector(y) ) )  stop("Arguments of wrong type")
   if ( (is.matrix(y) && nrow(y) != nrow(x)) || ( is.vector(y) && length(y) != nrow(x) )  )
     stop("Argument dimensions do not match")
   Sx <- var(x)
   Sxy <- var(x,y)
   Sy  <- var(y)
   if (is.vector(y) || ncol(y) == 1)  {
      Sy <- as.numeric(Sy)
      if (Sy < .Machine$double.eps) 
        stop("No associations can be found because the y vector is constant within the limits of machine precision")
      else H <- Sxy %*% t(Sxy)/Sy
      ry <- 1
   }
   else {  
      ry <- qr(Sy)$rank
      if (ry == ncol(y)) H <- Sxy %*% solve(Sy,t(Sxy))
      else  {
         if (ry == 1) warning("The y matrix has only one lineary independent column.\n")
         else warning("The y matrix has only ", ry, " lineary independent columns.\n")
         if (ry < ncol(x)) warning(" The expected rank of the H matrix was  adjusted accordingly.\n")
         library(MASS)
	 H <- Sxy %*% ginv(Sy) %*% t(Sxy)
      }
   }
   res <- list(mat=Sx,H=H,r=min(ncol(x),ry),call=match.call())
   res
} 

ldaHmat.default <- function(x,grouping,...)
{
  if  ( !is.matrix(x) || !is.factor(grouping) ) stop("Arguments of wrong type")
  n <- nrow(x)
  if (n != length(grouping)) stop("Argument dimensions do not match")
  nk <- table(grouping)
  k <- nrow(nk)
  T <- (n-1)*var(x)
  Si <- by(x,grouping,var)
  H <- T
  for (i in 1:k) H <- H - (nk[[i]]-1)*Si[[i]]
  res <- list(mat=T,H=H,r=min(ncol(x),k-1),call=match.call())
  res
}

glhHmat.default <- function(x,A,C,...)
{
   if  ( !is.matrix(x) || !is.matrix(A) || ( !is.matrix(C) && !is.vector(C) ) )   stop("Arguments of wrong type")
   if (is.vector(C))  C <- matrix(C,1,length(C))
   if (  nrow(A) != nrow(x) || ncol(C) != ncol(A) )  stop("Argument dimensions do not match")
   library(MASS)
   a0  <- A %*% ginv(t(A) %*% A) 
   Pa  <- a0 %*% t(A)
   ac  <- a0 %*% t(C)
   Pac <- ac %*% ginv(t(ac) %*% ac) %*% t(ac)
   T <- t(x) %*% (diag(1,nrow(x)) - Pa + Pac ) %*% x
   H   <- t(x) %*% Pac %*% x
   res <- list(mat=T,H=H,r=qr(C)$rank,call=match.call())
   res
}

lmHmat.data.frame <- function(x,y,...)
{
   res <- lmHmat.default(data.matrix(x),as.matrix(y))
   res$call <- match.call()
   res
}

ldaHmat.data.frame <- function(x,grouping,...)
{
   res <- ldaHmat.default(data.matrix(x),grouping)
   res$call <- match.call()
   res
}

glhHmat.data.frame <- function(x,A,C,...)
{
   if (is.vector(C)) res <- glhHmat.default(data.matrix(x),as.matrix(A),C)
   else res <- glhHmat.default(data.matrix(x),as.matrix(A),as.matrix(C))
   res$call <- match.call()
   res
}

lmHmat.formula <- function(formula,data=NULL,...)
{
   m <- match.call()
   if (is.matrix(eval(m$data,sys.parent()))) m$data <- as.data.frame(data)
   m[[1]] <- as.name("model.frame")
   m <- eval(m,sys.parent())
   Terms <- attr(m,"terms")
   y <- model.extract(m,"response")
   x <- model.matrix(Terms,m)
   xint <- match("(Intercept)",dimnames(x)[[2]],nomatch=0)
   if (xint>0) x <- x[,-xint,drop=F]
   res <- lmHmat.default(x,y)
   res$call <- match.call()
   res 
}

ldaHmat.formula <- function(formula,data=NULL,...)
{
   m <- match.call()
   if (is.matrix(eval(m$data,sys.parent()))) m$data <- as.data.frame(data)
   m[[1]] <- as.name("model.frame")
   m <- eval(m,sys.parent())
   Terms <- attr(m,"terms")
   grouping <- model.extract(m,"response")
   x <- model.matrix(Terms,m)
   xint <- match("(Intercept)",dimnames(x)[[2]],nomatch=0)
   if (xint>0) x <- x[,-xint,drop=F]
   res <- ldaHmat.default(x,grouping)
   res$call <- match.call()
   res
}

glhHmat.formula <- function(formula,C,data=NULL,...)
{
   m <- match.call()
   m$C <- NULL
   if (is.matrix(eval(m$data,sys.parent()))) m$data <- as.data.frame(data)
   m[[1]] <- as.name("model.frame")
   m <- eval(m,sys.parent())
   Terms <- attr(m,"terms")
   x <- model.extract(m,"response")
   A <- model.matrix(Terms,m)
   if (is.vector(C)) res <- glhHmat.default(as.matrix(x),A,C)
   else res <- glhHmat.default(as.matrix(x),A,as.matrix(C))
   res$call <- match.call()
   res
}


