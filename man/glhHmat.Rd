\name{glhHmat}
\alias{glhHmat}
\alias{glhHmat.default}
\alias{glhHmat.formula}
\alias{glhHmat.data.frame}
\title{Total and Effect Deviation Matrices for General Linear Hypothesis}
\description{Computes total an effect matrices of Sums of Squares and Cross-Product (SSCP) deviations for a general multivariate effect characterized by the violation of a linear hypothesis. These matrices may be used as input to the variable selection search routines \code{\link{anneal}}, \code{\link{genetic}} \code{\link{improve}} or \code{\link{leaps}}.}  

\usage{

\method{glhHmat}{default}(x,A,C,...)

\method{glhHmat}{data.frame}(x,A,C,...)

\method{glhHmat}{formula}(formula,C,data=NULL,...)

}
\arguments{
  \item{x}{A matrix or data frame containing the variables for which the SSCP matrix is to be computed.}
  \item{A}{A matrix or data frame containing a design matrix
specifying a linear model in which x is the response.} 
  \item{C}{A matrix or vector containing the coefficients of the
reference hypothesis.} 
  \item{formula}{A formula of the form \code{'x ~ A1 + A2 + ...'} That
is, the response is the set of variables whose subsets are to be
compared and the right hand side specifies the columns of the design
matrix.} 
  \item{data}{Data frame from which variables specified in 'formula'
are preferentially to be taken.}
  \item{...}{further arguments for the method.}  
}
\details{Consider a multivariate linear model \eqn{x = A Psi + U} and a
reference hypothesis \eqn{H0: C Psi = 0}, with \eqn{Psi} being a
matrix of unknown 
parameters and C a known coefficient matrix with rank \code{r}. It is well
known that, under classical Gaussian assumptions, \eqn{H_0} can be tested by
several increasing functions of the \code{r} positive 
eigenvalues of a product \eqn{T^{-1} H}, where \eqn{T} and \eqn{H} are
total and effect 
matrices of SSCP deviations associated with \eqn{H_0}. Furthermore, whether
or not the classical assumptions hold, the same eigenvalues can be
used to define descriptive indices that measure an "effect"
characterized by the violation of \eqn{H_0} (see reference [1] for further
details).  
Those SSCP matrices are given by \eqn{T = x'(I - P_{omega}) x} and \eqn{H =
x'(P_{Omega} - P_{omega}) x}, where I is an identity matrix and
\eqn{P_{Omega} = A(A'A)^-A'} , 

\deqn{P_{omega} = A(A'A)^-A' - A(A'A)^-C'[C(A'A)^-C']^-C(A'A)^-A'}

are projection matrices on the spaces spanned by the columns of A
(space \eqn{Omega}) and by the linear combinations of these columns that
satisfy the reference hypothesis (space \eqn{omega}). In these
formulae \eqn{M'} denotes the transpose of \eqn{M} and \eqn{M^-} a
generalized inverse. \code{glhHmat}  computes the \eqn{T} and \eqn{H}
matrices which then can be used as input to the 
search routines \code{\link{anneal}}, \code{\link{genetic}}
\code{\link{improve}} and \code{\link{leaps}} that try to select
subsets of x according to their contribution to the violation of \eqn{H_0}.}

  
\value{
A list with four items:
  \item{mat}{The total SSCP matrix}
  \item{H}{The effect SSCP matrix}
  \item{r}{The expected rank of the H matrix which equals the rank of
C. The true rank of H can be different from r if the x variables are
linearly dependent.} 
  \item{call}{The function call which generated the output.}
}

\references{[1] Duarte Silva. A.P. (2001). Efficient Variable
Screening for Multivariate Analysis, \emph{Journal of Multivariate
Analysis}, Vol. 76, 35-62.} 

\seealso{\code{\link{anneal}}, \code{\link{genetic}},
\code{\link{improve}}, \code{\link{leaps}}, \code{\link{lmHmat}},
\code{\link{ldaHmat}}.} 
\examples{
##----------------------------------------------------------------------------

##  The following examples create T and H matrices for different
## analysis of the MASS data set "crabs". This data 
##  records physical measurements on 200 specimens of Leptograpsus
## variegatus crabs observed on the shores of 
##  Western Australia. The crabs are classified by two factors, sex
## and sp (crab species as defined by its 
##  colour: blue or orange), with two levels each. The measurement
## variables include the carapace length (CL), c
##  the carapace width (CW), the size of the frontal lobe (FL) and the
## size of the rear width (RW). In the  
##  analysis provided, we assume that there is an interest in
## comparing the subsets of these variables measured in their original 
##  and logarithmic scales.        

library(MASS)
data(crabs)
lFL <- log(crabs$FL)
lRW <- log(crabs$RW)
lCL <- log(crabs$CL)
lCW <- log(crabs$CW)

# 1)  Create the T and H matrices associated with a linear
# discriminant analysis on the groups defined by the sp factor.  
# This call is equivalent to ldaHmat(sp ~ FL + RW + CL + CW  + lFL +
# lRW + lCL + lCW,crabs) 

glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sp,c(0,1),crabs)

# 2) Create the T and H matrices associated with a linear discriminant
# analysis on the groups defined by the sex factor. 
# This call is equivalent to ldaHmat(sex ~ FL + RW + CL + CW  + lFL +
# lRW + lCL + lCW,crabs) 

glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sex,c(0,1),crabs)

# 3)  Create the T and H matrices associated with a linear
# discriminant analysis on the groups 
# defined by all the combinations of the sp and sex factors

C <- matrix(0.,3,4)
C[row(C)+1 == col(C)] = 1.
glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sp*sex,C,crabs)

# 4)  Create the T and H matrices associated with an analysis
# of the interactions between the sp and sex factors

glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sp*sex,c(0,0,0,1),crabs)

# 5)  Create the T and H matrices associated with an analysis
# of the effect of the sp factor after controlling for sex

C <- matrix(0.,2,4)
C[1,3] = C[2,4] = 1.
glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sp*sex,C,crabs)

}
\keyword{manip}

