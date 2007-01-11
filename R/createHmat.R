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
   n <- nrow(x)
   p <- ncol(x)
   if (is.vector(y)) q <- 1
   else q <- ncol(y)
   if ( (is.matrix(y) && nrow(y) != n) || ( is.vector(y) && length(y) != n )  )
     stop("Argument dimensions do not match")

   Sx <- var(x)
   dx <- x - matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   if (q==1) {
	dy <- y - rep(mean(y),n)
	uy <- dy / sqrt(sum(dy^2))
   }	
   else  {
	dy <- y - matrix(rep(apply(y,2,mean),n),n,q,byrow=T)
   	uy <- svd(dy,nv=0)$u
   }
   h <- t(uy) %*% dx
   H <- t(h) %*% h / (n-1)

  res <- list(mat=Sx,H=H,r=min(p,q),call=match.call())
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

   r <- qr(C)$rank
   rA <- qr(A)$rank
   if ( rA <= r) stop("There are not enough linearly independent columns in the desing matrix (A)")
   svdA <- svd(A)
   a <-  t(svdA$u) %*% x 
   E <- t(x) %*% x - t(a) %*% a
   M <- svdA$u[,1:rA] %*% diag(svdA$d[1:rA]^-1) %*% t(svdA$v[,1:rA]) %*% t(C)
   h <- t(svd(M,nv=0)$u) %*% x
   H <- t(h) %*% h

   res <- list(mat=H+E,H=H,r=r,call=match.call())
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

