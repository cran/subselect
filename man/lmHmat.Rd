\name{lmHmat}
\alias{lmHmat}
\alias{lmHmat.default}
\alias{lmHmat.formula}
\alias{lmHmat.data.frame}
\title{Total and Effect Deviation Matrices for Linear Regression and Canonical Correlation Analysis}
\description{Computes total an effect matrices of Sums of Squares and Cross-Product (SSCP) deviations, divided by a normalizing constant, in linear regression or canonical correlation analysis. These matrices may be used as input to the variable selection search routines \code{\link{anneal}}, \code{\link{genetic}} \code{\link{improve}} or \code{\link{leaps}}.}  

\usage{

\method{lmHmat}{default}(x,y,...)

\method{lmHmat}{data.frame}(x,y,...)

\method{lmHmat}{formula}(formula,data=NULL,...)

}
\arguments{
  \item{x}{A matrix or data frame containing the variables for which
the SSCP matrix is to be computed.} 
  \item{y}{A matrix or data frame containing the set of fixed
variables, the association of x is to be measured with.} 
  \item{formula}{A formula of the form \code{'y ~ x1 + x2 +
...'}. That is, the response is the set of fixed variables and the
right hand side specifies the variables whose subsets are to be
compared.} 
  \item{data}{Data frame from which variables specified in 'formula'
are preferentially to be taken.} 
  \item{...}{further arguments for the method.} 
}
\details{Let x and y be two different groups of linearly independent
variables observed on the same set of data units. It is well known
that the association between x and y can be measured by their squared
canonical correlations which may be found as the positive eigenvalues
of certain matrix products. In particular, if \eqn{T_x} and
\eqn{H_{x/y}} denote SSCP 
matrices of deviations from the mean, respectively for the original x
variables (\eqn{T_x}) and for their orthogonal projections onto the space
spanned by the y's (\eqn{H_{x/y}}), then the positive eigenvalues of
\eqn{T_x^{-1}H_{x/y}} equal the squared correlations between x and
y. Alternatively 
these correlations could also be found from \eqn{T_y^{-1} H_{y/x}} but here,
assuming a goal of comparing x's subsets for a given fixed set of y's,
we will focus on the former product. \code{lmHmat} computes a scaled version
of \eqn{T_x} and \eqn{H_{x/y}} such that \eqn{T_x} is converted into a
covariance matrix. These matrices can be used as input to the search routines
\code{\link{anneal}}, \code{\link{genetic}} \code{\link{improve}} and
\code{\link{leaps}} that try to select x subsets based on several
functions of their squared correlations with y. We note that when
there is only one variable in the y set, this is equivalent to
selecting predictors for linear regression based on the traditional
coefficient of determination. }    
  
\value{
A list with four items:
  \item{mat}{The total SSCP matrix divided by nrow(x)-1}
  \item{H}{The effect SSCP matrix divided by nrow(x)-1}
  \item{r}{The expected rank of the H matrix which, under the
assumption of linear independence, equals the minimum between the
number of variables in the x and y sets. The true rank of \code{H} can be
different from \code{r} if the linear independence condition fails.} 
  \item{call}{The function call which generated the output.}
}
\seealso{\code{\link{anneal}}, \code{\link{genetic}},
\code{\link{improve}}, \code{\link{leaps}}, \code{\link{lm}}.} 
\examples{
##------------------------------------------------------------------

## 1)  An example of subset selection in the context of Multiple
## Linear Regression. Variable 5 (average price) in the Cars93 MASS
## library is to be regressed on 13 other  variables.  The goal is to
## compare subsets of these 13 variables according to their  ability
## to predict car prices. 

library(MASS)
data(Cars93)
lmHmat(Cars93[c(7:8,12:15,17:22,25)],Cars93[5])

## 2)  An example of subset selection in the context of Canonical
## Correlation Analysis. Two groups of variables within the Cars93
## MASS library data set are compared. The first group (variables 4th,
## 5th and 6th) relates to price, while the second group is formed by 13
## variables that describe several technical car specifications. The
## goal is to select subsets of the second group that are optimal in
## terms of preserving the canonical correlations with the variables in
## the first group (Warning: the 3-variable "response" group is kept
## intact; subset selection is to be performed only in the 13-variable
## group).  

library(MASS)
data(Cars93)
lmHmat(Cars93[c(7:8,12:15,17:22,25)],Cars93[4:6])
}
\keyword{manip}


