\name{anneal}
\alias{anneal}
\title{SIMULATED ANNEALING SEARCH FOR AN OPTIMAL k-VARIABLE SUBSET} 
\description{Given a set of variables,  a Simulated Annealing
algorithm seeks a k-variable subset which is
optimal, as a surrogate for the whole set, with respect to a given criterion. 
}
\details{An initial k-variable subset (for k ranging from \emph{kmin} to \emph{kmax}) 
of a full set of p (p not exceeding 300) variables is randomly
selected and passed on to a Simulated Annealing algorithm. 
The algorithm then selects a random subset in the neighbourhood of the
current subset (neighbourhood of a subset S being defined as the
family of all k-variable subsets which differ from S by a
single variable), and decides whether to replace the current subset 
according to the Simulated Annealing rule, i.e., either (i) always,
if the alternative subset's value of the criterion is higher; or (ii) with
probability \eqn{$\exp^{\frac{ac-cc}{t}}$}{exp((ac-cc)/t)} 
if the alternative subset's value of the
criterion (ac) is lower than that of the current solution (cc), where
the parameter t (temperature) decreases throughout the
iterations of the algorithm. For each
cardinality k, the stopping criterion for the
algorithm is the number of iterations (\emph{niter}) which is controlled by the
user. Also controlled by the user are the initial temperature (\emph{temp}) and the
rate of geometric cooling of the temperature (\emph{cooling}). 

Optionally, the best k-variable subset produced by Simulated Annealing may be
passed as input to a restricted local search algorithm, for possible further
improvement.

The user may force variables to be included and/or excluded from the k-subsets.

For each cardinality k, the total number of calls to the procedure which computes the criterion
values is \emph{nsol} x (\emph{niter} + 1). These calls are the dominant computational
effort in each iteration of the algorithm.  

In order to improve computation times, the bulk of computations are
carried out by a Fortran routine which is called by this
function. Further details about the Simulated Annealing algorithm can
be found in Reference 1 and in the comments to the Fortran code (in
the \code{src} subdirectory for this package). For p > 300, it is 
necessary to change the declarative statements in the Fortran code. 
}
\usage{anneal( mat, kmin, kmax = kmin, nsol = 1, niter = 1000, exclude
= NULL, include = NULL, improvement = TRUE, printfile = FALSE, iseed =
c(sample(0:4095,3),sample(0:2047,1)*2+1), cooling = 0.05,  temp = 1,
criterion = "RM", pcindices = 1:kmax))} 
\arguments{
  \item{mat}{a covariance or correlation matrix of the variables from
  which the k-subset is to be selected.}
  \item{kmin}{the cardinality of the smallest subset that is wanted.}
  \item{kmax}{the cardinality of the largest subset that is wanted.}
  \item{nsol}{the number of initial subsets (runs of the algorithm).}
  \item{niter}{the number of iterations of the algorithm for each
  initial subset.}  
  \item{exclude}{a vector of variables (referenced by their row/column
  numbers in matrix \emph{mat}) that are to be forcibly excluded from the subsets.}
  \item{include}{a vector of variables (referenced by their row/column
  numbers in matrix \emph{mat}) that are to be forcibly included in the subsets.}
  \item{improvement}{a logical variable indicating whether or not the
  best final subset (for each cardinality) is to be passed as input to a
  local improvement algorithm (see \code{\link{improve}}).}
  \item{printfile}{a logical variable indicating whether results in
  print-friendly format are to be written to a file.}
  \item{iseed}{a vector of 4 integers (as requested by the LAPACK
  routine SLARUV) for the seeds of the random number generation.}
  \item{cooling}{variable in the ]0,1[ interval indicating the rate of
  geometric cooling for the Simulated Annealing algorithm.}
  \item{temp}{positive variable indicating the initial
  temperature for the Simulated Annealing algorithm.}
  \item{criterion}{Indicates which criterion is to be used in judging
  the quality of the subsets. Currently, only the RM, RV and GCD criteria
  are supported (see References, \code{\link{rm.coef}},
  \code{\link{rv.coef}} and \code{\link{gcd.coef}} for further details).}
  \item{pcindices}{a vector of ranks of Principal Components that are to be
  used for comparison with the k-variable subsets (for the GCD
  criterion only, see \code{\link{gcd.coef}}).}
}
\value{A list with four items:
   \item{subsets}{An \emph{nsol} x \emph{kmax} x length(\emph{kmin}:\emph{kmax}) 3-dimensional
   array, giving for each cardinality (dimension 3) and each solution
   (dimension 1) the list of variables (referenced by their row/column
  numbers in matrix \emph{mat}) in the subset (dimension 2). (For
  cardinalities  smaller than \emph{kmax}, the extra positions are set to zero).}
   \item{values}{An \emph{nsol} x length(\emph{kmin}:\emph{kmax}) matrix, giving for each
   cardinality (columns), the criterion values of the \emph{nsol} (rows) subsets obtained.}
   \item{bestvalues}{A length(\emph{kmin}:\emph{kmax}) vector giving the best values
   of the criterion obtained for each cardinality. If \emph{improvement} is
   TRUE, these values result from the final restricted local
   search algorithm (and may therefore exceed the largest value for
   that cardinality in \emph{values}).}
   \item{bestsets}{A length(\emph{kmin}:\emph{kmax}) x \emph{kmax} matrix, giving, for each
   cardinality (rows), the variables (referenced by their row/column numbers in matrix \emph{mat}) in the
   best k-subset that was found.}
}
\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}}, \code{\link{genetic}}, \code{\link{improve}}.}
\references{
1) Cadima, J., Cerdeira, J. Orestes and Minhoto, M. A Comparative
Study of Algorithms for Variable Selection in Multivariate
Statistics, submitted for publication in \emph{American Journal of
Mathematical and Management Sciences}.

2) Cadima, J. and Jolliffe, I.T. (2001). Variable Selection and the
Interpretation of Principal Subspaces, \emph{Journal of Agricultural,
Biological and Environmental Statistics}, Vol. 6, 62-79.
}

\examples{

# For illustration of use, a small data set with very few iterations
# of the algorithm.  (iseed is used to ensure the same solution,
# given the random nature of the algorithm).

data(swiss)
anneal(cor(swiss),2,3,nsol=4,niter=10,criterion="RM",iseed=c(2345,3276,1568,3187))
## $subsets
## , , Card.2
## 
##            Var.1 Var.2 Var.3
## Solution 1     3     6     0
## Solution 2     1     2     0
## Solution 3     3     6     0
## Solution 4     3     6     0
##
## , , Card.3
##
##            Var.1 Var.2 Var.3
## Solution 1     3     5     6
## Solution 2     1     2     5
## Solution 3     1     2     5
## Solution 4     4     5     6
##
##
## $values
##               card.2    card.3
## Solution 1 0.8016409 0.8769672
## Solution 2 0.7945390 0.8791856
## Solution 3 0.8016409 0.8791856
## Solution 4 0.8016409 0.9043760
##
## $bestvalues
##    Card.2    Card.3 
## 0.8016409 0.9043760 
##
## $bestsets
##        Var.1 Var.2 Var.3
## Card.2     3     6     0
## Card.3     4     5     6
##


#
#
# Excluding variable number 6 from the subsets.
# 

data(swiss)
anneal(cor(swiss),2,3,nsol=4,niter=10,criterion="RM",exclude=c(6),iseed=c(2345,3276,1568,3187))
## $subsets
## , , Card.2
##
##            Var.1 Var.2 Var.3
## Solution 1     4     5     0
## Solution 2     4     5     0
## Solution 3     4     5     0
## Solution 4     4     5     0
##
## , , Card.3
##
##            Var.1 Var.2 Var.3
## Solution 1     1     2     5
## Solution 2     1     2     5
## Solution 3     1     2     5
## Solution 4     1     2     5
##
##
## $values
##               card.2    card.3
## Solution 1 0.7982296 0.8791856
## Solution 2 0.7982296 0.8791856
## Solution 3 0.7982296 0.8791856
## Solution 4 0.7982296 0.8791856
##
## $bestvalues
##    Card.2    Card.3 
## 0.7982296 0.8791856 
##
## $bestsets
##        Var.1 Var.2 Var.3
## Card.2     4     5     0
## Card.3     1     2     5
}
\keyword{manip}
