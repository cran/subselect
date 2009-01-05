\name{zeta2.coef}

\alias{zeta2.coef}

\title{Computes the Zeta squared coefficient for variable subset selection} 
\description{
Computes the Zeta squared  index of "effect magnitude".  
The maximization of this criterion is equivalent to the maximization of the traditional test statistic, the Lawley-Hotelling  trace.

}

\usage{zeta2.coef(mat, H=NULL, r=0, indices,
tolval=10*.Machine$double.eps, tolsym=1000*.Machine$double.eps)} 

\arguments{
  \item{mat}{the Variance or Total sums of squares and products matrix for the full data set.}

  \item{H}{the Effect description sums of squares and products matrix (defined in the same way as the \code{mat} matrix).}

  \item{r}{the Expected rank of the H matrix. See the \code{Details} below.}

  \item{indices}{a numerical vector, matrix or 3-d array of integers
  giving the indices of the variables in the subset. If a matrix is
  specified, each row is taken to represent a different
  \emph{k}-variable subset. If a 3-d array is given, it is assumed
  that the third dimension corresponds to different cardinalities.}

  \item{tolval}{the tolerance level to be used in checks for 
ill-conditioning and positive-definiteness of the 'total' and
'effects' (H) matrices. Values smaller than \code{tolval} are
considered equivalent to zero.}

 \item{tolsym}{the tolerance level for symmetry of the
   covariance/correlation/total matrix and for the effects (\code{H})
   matrix. If corresponding matrix entries differ by more than this
   value, the input matrices will be considered asymmetric and execution
   will be aborted. If corresponding entries are different, but by less
   than this value, the input matrix will be replaced by its symmetric
   part, i.e., input matrix A becomes (A+t(A))/2.}
}

\details{
Different kinds of statistical methodologies are considered within the
framework, of a multivariate linear model:
\deqn{X = A \Psi + U}{X = A B + U} where \eqn{X} is the (nxp) data
matrix of original variables, \eqn{A} is a known (nxp) design matrix, 
\eqn{\Psi}{B} an (qxp) matrix of unknown parameters and \eqn{U} an
(nxp) matrix of residual vectors. The \eqn{\zeta^2}{zeta2}  index is related to the traditional test statistic (Lawley-Hotelling  trace) and 
measures the contribution of each subset to an Effect characterized by
the violation of a linear hypothesis of the form 
\eqn{C \Psi = 0}{C B = 0}, where \eqn{C} is a known cofficient matrix
of rank r. The Lawley-Hotelling  trace is given by: 
\eqn{V=tr(HE^{-1})} where \eqn{H} is the Effect matrix and \eqn{E} is
the Error matrix. 
 The index \eqn{\zeta^2}{zeta2} is related to Lawley-Hotelling  trace
(\eqn{V}) by:  \deqn{\zeta^2 =\frac{V}{V+r}}{zeta^2 = V/(V+r)} where
\eqn{r} is the rank of \eqn{H} matrix.




The fact that \code{indices} can be a matrix or 3-d array allows for
the computation of the \eqn{\zeta^2} values of subsets produced by the search
functions \code{\link{anneal}}, \code{\link{genetic}}, 
\code{\link{improve}} and
\code{\link{leaps}} (whose output option \code{$subsets} are
matrices or 3-d arrays), using a different criterion (see the example
below).
}

\value{
The value of the \eqn{\zeta^2}{zeta2} coefficient.
}

\references{
 Duarte Silva, A.P. (2001) Efficient Variable Screening for Multivariate Analysis,
\emph{Journal of Multivariate Analysis}, Vol. 76, 35-62.

}

\examples{

## ---------------------------------------------------------------

## 1) A Linear Discriminant Analysis example with a very small data set. 
## We considered the Iris data and three groups, 
## defined by species (setosa, versicolor and virginica).
## Both the 'iris' and the 'iris3' R objects are used, for convenience in  
## computing the Total matrix (T=irismat) and the Effects matrix (H=irish).

data(iris)
data(iris3) 
irismat <-149*var(iris[,1:4])
weights<-c(49,49,49)
irish <- irismat - matrix(weights\%*\% t(apply(iris3,3,var)),4,4)
zeta2.coef(irismat,H=irish,r=2,c(1,3))
## [1] 0.9211501

## ---------------------------------------------------------------

## 2) An example computing the value of the zeta_2 criterion for two  
## subsets produced when the anneal function attempted to optimize  
## the ccr1_2 criterion (using an absurdly small number of iterations).

ccr1results<-anneal(irismat,2,nsol=2,niter=2,criterion="ccr12",H=irish,r=2)
zeta2.coef(irismat,H=irish,r=2,ccr1results$subsets)

##              Card.2
##Solution 1 0.9105021
##Solution 2 0.9161813

## ---------------------------------------------------------------

}

\keyword{manip}
