\name{genetic}
\alias{genetic}
\title{GENETIC ALGORITHM FOR AN OPTIMAL k-VARIABLE SUBSET} 
\description{Given a set of variables,  a Genetic Algorithm
algorithm seeks a k-variable subset which is
optimal, as a surrogate for the whole set, with respect to a given criterion. 
}
\details{
For each cardinality k (with k ranging from \emph{kmin} to \emph{kmax}),
an initial population of \emph{popsize} k-variable subsets is randomly
selected from a full set of p (p not exceeding 300) variables. 
In each iteration, \emph{popsize}/2 couples
are formed from among the population and each couple generates a child
(a new k-variable subset)
which inherits properties of its parents (specifically, it inherits
all variables common to both parents and a random selection of
variables in the symmetric difference of its parents' genetic makeup).
Each offspring may optionally undergo a mutation (in the form of a
local improvement algorithm (see \code{\link{improve}}). The parents
and offspring are ranked according to their criterion value, and the
best \emph{popsize} of these k-subsets will make up the next
generation, which is used as the current population in the subsequent iteration.

The stopping rule for the algorithm is the number of generations (\emph{nger}).

Optionally, the best k-variable subset produced by the Genetic Algorithm may be
passed as input to a restricted local search algorithm, for possible further
improvement.

The user may force variables to be included and/or excluded from the k-subsets.

For each cardinality k, the total number of calls to the procedure which computes the criterion
values is \eqn{popsize + nger} x \eqn{popsize/2}. These calls are the dominant computational
effort in each iteration of the algorithm.  

In order to improve computation times, the bulk of computations are
carried out by a Fortran routine which is called by this
function. Further details about the Genetic Algorithm can
be found in Reference 1 and in the comments to the Fortran code (in
the \code{src} subdirectory for this package). For p > 300, it is 
necessary to change the declarative statements in the Fortran code. 
}
\usage{genetic( mat, kmin, kmax = kmin, popsize = 100, nger = 20,
mutate = FALSE, maxclone = 5, exclude = NULL, include = NULL,
improvement = TRUE, printfile = FALSE,
iseed = c(sample(0:4095,3),sample(0:2047,1)*2+1), criterion = "RM",
pcindices = 1:kmax)
}
\arguments{
  \item{mat}{a covariance or correlation matrix of the variables from
  which the k-subset is to be selected.}
  \item{kmin}{the cardinality of the smallest subset that is wanted.}
  \item{kmax}{the cardinality of the largest subset that is wanted.}
  \item{popsize}{integer variable indicating the size of the population.}
  \item{nger}{integer variable giving the number of generations for
  which the genetic algorithm will run.}  
  \item{mutate}{logical variable indicating whether each  child
  undergoes a mutation. By default, FALSE.}
  \item{maxclone}{integer variable specifying the maximum number of
  identical replicates (clones) of individuals that is acceptable in
  the population. Serves to ensure that the population has sufficient
  genetic diversity, which is necessary to enable the algorithm to
  complete the specified number of generations.}
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
  \item{criterion}{Indicates which criterion is to be used in judging
  the quality of the subsets. Currently, only the RM, RV and GCD criteria
  are supported (see References, \code{\link{rm.coef}},
  \code{\link{rv.coef}} and \code{\link{gcd.coef}} for further details).}
  \item{pcindices}{a vector of ranks of Principal Components that are to be
  used for comparison with the variable subsets (for the GCD
  criterion only, see \code{\link{gcd.coef}}).}
}
\value{A list with four items:
   \item{subsets}{A \emph{popsize} x \emph{kmax} x length(\emph{kmin}:\emph{kmax}) 3-dimensional
   array, giving for each cardinality (dimension 3) and each subset in
   the final population  (dimension 1) the list of variables (referenced by their row/column
  numbers in matrix \emph{mat}) in the subset (dimension 2). (For
  cardinalities  smaller than \emph{kmax}, the extra positions are set to zero).}
   \item{values}{A \emph{popsize} x length(\emph{kmin}:\emph{kmax}) matrix, giving for each
   cardinality (columns), the (ordered) criterion values of the
   \emph{popsize} (rows) subsets in the final generation.}
   \item{bestvalues}{A length(\emph{kmin}:\emph{kmax}) vector giving the best values
   of the criterion obtained for each cardinality. If \emph{improvement} is
   TRUE, these values result from the final restricted local
   search algorithm (and may therefore exceed the largest value for
   that cardinality in \emph{values}).}
   \item{bestsets}{A length(\emph{kmin}:\emph{kmax}) x \emph{kmax} matrix, giving, for each
   cardinality (rows), the variables (referenced by their row/column
  numbers in matrix \emph{mat}) in the
   best k-subset that was found.}
}
\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}}, \code{\link{anneal}}, \code{\link{improve}}.}
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

data(swiss)
genetic(cor(swiss),3,4,popsize=10,nger=5,iseed=c(2345,3276,1568,3187),criterion="Rv")
  
##  For cardinality k= 4
##  Generation number 4 does not
##  have sufficient genetic diversity
##  for acceptable levels of consanguinity
##  (couples differing by at least 2 genes).
##  
##  This situation may be overcome by 
##  reducing the maximum acceptable number 
##  of clones (maxclone).
##  
##  
##  Best criterion value found so far:  0.955714479
##  
## $subsets
##             Var.1 Var.2 Var.3
## Solution 1      1     2     3
## Solution 2      1     2     5
## Solution 3      3     4     6
## Solution 4      1     4     5
## Solution 5      3     4     5
## Solution 6      3     4     5
## Solution 7      3     4     5
## Solution 8      2     4     5
## Solution 9      1     3     4
## Solution 10     1     3     4
##
## $values
##  Solution 1  Solution 2  Solution 3  Solution 4  Solution 5  Solution 6 
##  0.9141995   0.9098502   0.9034868   0.9023910   0.9020271   0.9020271 
##  Solution 7  Solution 8  Solution 9 Solution 10 
##  0.9020271   0.8982510   0.8940945   0.8940945 
##
## $bestvalues
##    Card.3 
## 0.9141995 
##
## $bestsets
## Var.1 Var.2 Var.3 
##     1     2     3 
##
}
\keyword{manip}
