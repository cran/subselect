\name{full.k.search}
\alias{full.k.search}
\title{VALUES OF THE RM, GCD AND RV COEFFICIENTS FOR ALL
k-VARIABLE SUBSETS OF A DATA SET}
\description{
Computes the values of the GCD, RV and RM
coefficients for all \emph{k}-variable subsets of a given data
set. Outputs the optimal values and subsets, for each of those
coefficients, and, optionally, the full results.
}
\usage{full.k.search(mat, k, print.all = FALSE, file="")}
\arguments{
  \item{mat}{the full data set's covariance (or correlation) matrix.}
  \item{k}{the cardinality of the variable subsets that are wanted.}
  \item{print.all}{if TRUE, prints out the values of the three
coefficients for all k-variable subsets.}
  \item{file}{the file where the complete results for all
k-variable subsets will be written, if \code{print.all}=TRUE.}
}
\details{
Generates all \emph{k}-variable subsets of the \emph{p}-variable data
set defined by 
\code{mat}. For each subset, computes the values of the GCD, RM and RV
coefficients (comparing with the first \emph{k} Principal Components
when computing \code{\link{gcd.coef}}). When \code{print.all}=FALSE, only
the optimal values and subsets, for each coefficient, are produced in
standard output. If \code{print.all}=TRUE, the full results are written in
the file specified in \code{file}. 
}
\section{Warning}{This function is unusable even for moderatly
large data sets. For such data sets, consider using the SSCMA software, written and
presented by Pedro Duarte Silva in \emph{Discarding Variables in Principal
Component Analysis: Algorithms for all-subset comparisons},
WP-00-002, July 2000, Universidade Católica Portuguesa, \email{psilva@porto.ucp.pt}.}
\value{
A list with the following items:
\item{rmmax}{The maximum value of the RM coefficient;}
\item{rmind}{The indices of the \emph{k} variables in the optimal
subset for the RM coefficient;}  
\item{gcdmax}{The maximum value of the GCD coefficient;} 
\item{gcdind}{The indices of the \emph{k} variables in the optimal
subset for the GCD coefficient;} 
\item{rvmax}{The maximum value of the RV coefficient;} 
\item{rvind}{The indices of the \emph{k} variables in the optimal
subset for the RV coefficient.} 

If \code{print.all}=TRUE, the full results for all \emph{k}-variable
subsets will be written to \code{file}, one subset per line.
}
\seealso{\code{\link{gcd.coef}},\code{\link{rm.coef}},\code{\link{rv.coef}}} 
\references{
 Cadima, J. and Jolliffe, I.T. (2001), "Variable Selection and the
Interpretation of Principal Subspaces", \emph{Journal of Agricultural,
Biological and Environmental Statistics}, Vol. 6, 62-79.
}
\examples{
data(iris3) 
x<-iris3[,,1]
full.k.search(cor(x),k=2)
## $rmmax
## [1] 0.8233218
## 
## $rmindices
## [1] 2 3
##
## $gcdmax
## [1] 0.7821987
##
## $gcdindices
## [1] 2 3
##
## $rvmax
## [1] 0.8453146
##
## $rvindices
## [1] 1 4
}
\keyword{manip}
