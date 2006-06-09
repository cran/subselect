\name{ldaHmat}
\alias{ldaHmat}
\alias{ldaHmat.default}
\alias{ldaHmat.formula}
\alias{ldaHmat.data.frame}
\title{Total and Between-Group Deviation Matrices in Linear Discriminant Analysis}
\description{Computes total and between-group matrices of Sums of Squares and Cross-Product (SSCP) deviations in linear discriminant analysis. These matrices may be used as input to the variable selection search routines \code{\link{anneal}}, \code{\link{genetic}} \code{\link{improve}} or \code{\link{leaps}}.}  

\usage{

\method{ldaHmat}{default}(x,grouping,...)

\method{ldaHmat}{data.frame}(x,grouping,...)

\method{ldaHmat}{formula}(formula,data=NULL,...)

}
\arguments{
  \item{x}{A matrix or data frame containing the discriminators for
which the SSCP matrix is to be computed.} 
  \item{grouping}{A factor specifying the class for each observation.}
  \item{formula}{A formula of the form \code{'groups ~ x1 + x2 + ...'}
That is, the response is the grouping factor and the right hand side
specifies the (non-factor) discriminators.} 
  \item{data}{Data frame from which variables specified in 'formula'
are preferentially to be taken.}  
  \item{...}{further arguments for the method.}
}
  
\value{
A list with four items:
  \item{mat}{The total SSCP matrix}
  \item{H}{The between-groups SSCP matrix}
  \item{r}{The expected rank of the H matrix which equals the minimum
between the number of discriminators and the number of groups minus
one. The true rank of H can be different from r if the discriminators
are linearly dependent.} 
  \item{call}{The function call which generated the output.}
}
\seealso{\code{\link{anneal}}, \code{\link{genetic}},
\code{\link{improve}}, \code{\link{leaps}}, \code{\link{lda}}.} 
\examples{
##--------------------------------------------------------------------

## An example with a very small data set. We consider the Iris data
## and three groups, defined by species (setosa, versicolor and
## virginica). 

data(iris)
ldaHmat(iris[1:4],iris$Species)

}
\keyword{manip}

