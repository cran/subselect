\name{genetic}

\alias{genetic}

\title{Genetic Algorithm searching for an optimal k-variable subset} 

\description{Given a set of variables,  a Genetic Algorithm
algorithm seeks a k-variable subset which is
optimal, as a surrogate for the whole set, with respect to a given criterion. 
}

\details{
For each cardinality k (with k ranging from \code{kmin} to \code{kmax}),
an initial population of \code{popsize} k-variable subsets is randomly
selected from a full set of p (p not exceeding 300) variables. 
In each iteration, \code{popsize}/2 couples
are formed from among the population and each couple generates a child
(a new k-variable subset)
which inherits properties of its parents (specifically, it inherits
all variables common to both parents and a random selection of
variables in the symmetric difference of its parents' genetic makeup).
Each offspring may optionally undergo a mutation (in the form of a
local improvement algorithm -- see function \code{\link{improve}}),
with a user-specified probability. The parents
and offspring are ranked according to their criterion value, and the
best \code{popsize} of these k-subsets will make up the next
generation, which is used as the current population in the subsequent
iteration. 

The stopping rule for the algorithm is the number of generations (\code{nger}).

Optionally, the best \emph{k}-variable subset produced by the Genetic
Algorithm may be passed as input to a restricted local improvement
algorithm, for possible further improvement (see function
\code{\link{improve}}). 

The user may force variables to be included and/or excluded from the
\emph{k}-subsets, and may specify an initial population.

For each cardinality \emph{k}, the total number of calls to the
procedure which computes the criterion 
values is \eqn{popsize + nger} x \eqn{popsize/2}. These calls are the
dominant computational effort in each iteration of the algorithm.  

In order to improve computation times, the bulk of computations are
carried out by a Fortran routine. Further details about the Genetic
Algorithm can 
be found in Reference 1 and in the comments to the Fortran code (in
the \code{src} subdirectory for this package). For p > 300, it is 
necessary to change the declarative statements in the Fortran code. 
}

\usage{genetic( mat, kmin, kmax = kmin, popsize = 100, nger = 100,
mutate = FALSE, mutprob = 0.01, maxclone = 5, exclude = NULL,
include = NULL, improvement = TRUE, setseed= FALSE, criterion = "RM",
pcindices = 1:kmax, initialpop=NULL )
}

\arguments{
  \item{mat}{a covariance or correlation matrix of the variables from
  which the k-subset is to be selected.}

  \item{kmin}{the cardinality of the smallest subset that is wanted.}

  \item{kmax}{the cardinality of the largest subset that is wanted.}

  \item{popsize}{integer variable indicating the size of the
  population.}

  \item{nger}{integer variable giving the number of generations for
  which the genetic algorithm will run.}  

  \item{mutate}{logical variable indicating whether each  child
  undergoes a mutation, with probability \code{mutprob}. By default, FALSE.}

  \item{mutprob}{variable giving the probability of each  child
  undergoing a mutation, if \code{mutate} is TRUE. By default, 0.01.
  High values slow down the algorithm considerably and tend to
  replicate the same solution.}

  \item{maxclone}{integer variable specifying the maximum number of
  identical replicates (clones) of individuals that is acceptable in
  the population. Serves to ensure that the population has sufficient
  genetic diversity, which is necessary to enable the algorithm to
  complete the specified number of generations. However, even maxclone=0
  does not guarantee that there are no repetitions: only the offspring 
  of couples are tested for clones. If any such clones are rejected, they  
  are replaced by a k-variable subset chosen at random, without any
  further clone tests.}

  \item{exclude}{a vector of variables (referenced by their row/column
  numbers in matrix \code{mat}) that are to be forcibly excluded from
  the subsets.} 

  \item{include}{a vector of variables (referenced by their row/column
  numbers in matrix \code{mat}) that are to be forcibly included in
  the subsets.} 

  \item{improvement}{a logical variable indicating whether or not the
  best final subset (for each cardinality) is to be passed as input to a
  local improvement algorithm (see function \code{\link{improve}}).}

  \item{setseed}{logical variable indicating whether to fix an initial 
  seed for the random number generator, which will be re-used in future
  calls to this function whenever setseed is again set to TRUE.}

  \item{criterion}{Character variable, which indicates which criterion
  is to be used in judging the quality of the subsets. Currently, only
  the RM, RV and GCD criteria are supported, and referenced as "RM",
  "RV" or "GCD" (see References, \code{\link{rm.coef}}, 
  \code{\link{rv.coef}} and \code{\link{gcd.coef}} for further details).}

  \item{pcindices}{a vector of ranks of Principal Components that are to be
  used for comparison with the variable subsets (for the GCD
  criterion only, see \code{\link{gcd.coef}}).}

  \item{initialpop}{vector, matrix or 3-d array of initial population
  for the genetic algorithm. If a \emph{single cardinality} is
  required, \code{initialpop} may be a \code{popsize} x \emph{k}
  matrix or a \code{popsize} x \emph{k} x 1 array (as produced by the
  \code{$subsets} output value of any of the 
  algorithm functions \code{anneal}, \code{genetic}, or
  \code{improve}). If \emph{more 
  than one cardinality} is requested, \code{initialpop} must be a
  \code{popsize x kmax x length(kmin:kmax)} 3-d array (as produced by the
  \code{$subsets} output value).

  If the \code{exclude} and/or \code{include} options are used,
  \code{initialpop} must also respect those requirements. } 
}

\value{A list with four items:

   \item{subsets}{A \code{popsize} x \code{kmax} x
   length(\code{kmin}:\code{kmax}) 3-dimensional array, giving for
   each cardinality (dimension 3) and each subset in the final
   population  (dimension 1) the list of variables (referenced by
   their row/column numbers in matrix \code{mat}) in the subset
   (dimension 2). (For cardinalities  smaller than \code{kmax}, the
   extra final positions are set to zero).} 

   \item{values}{A \code{popsize} x length(\code{kmin}:\code{kmax})
   matrix, giving for each cardinality (columns), the (ordered)
   criterion values of the \code{popsize} (rows) subsets in the final
   generation.} 

   \item{bestvalues}{A length(\code{kmin}:\code{kmax}) vector giving
   the best values of the criterion obtained for each cardinality. If
   \code{improvement} is TRUE, these values result from the final
   restricted local search algorithm (and may therefore exceed the
   largest value for that cardinality in \code{values}).}

   \item{bestsets}{A length(\code{kmin}:\code{kmax}) x \code{kmax}
   matrix, giving, for each cardinality (rows), the variables
   (referenced by their row/column numbers in matrix \code{mat}) in the
   best k-subset that was found.}
}

\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}}, \code{\link{anneal}}, \code{\link{improve}}.}

\references{
1) Cadima, J., Cerdeira, J. Orestes and Minhoto, M. A Comparative
Study of Algorithms for Variable Selection in Multivariate
Statistics, submitted for publication.

2) Cadima, J. and Jolliffe, I.T. (2001). Variable Selection and the
Interpretation of Principal Subspaces, \emph{Journal of Agricultural,
Biological and Environmental Statistics}, Vol. 6, 62-79.
}

\examples{
# For illustration of use, a small data set with very few iterations
# of the algorithm.  

data(swiss)
genetic(cor(swiss),3,4,popsize=10,nger=5,criterion="Rv")

## For cardinality k=
##[1] 4
## there is not enough genetic diversity in generation number 
##[1] 5
## for acceptable levels of consanguinity (couples differing by at
## least 2 genes). 
## [1]
## Try reducing the maximum acceptable number  of clones (maxclone) or
## increasing the population size (popsize) 
## [1]
## Best criterion value found so far:
##[1] 0.9590526
##$subsets
##            Var.1 Var.2 Var.3
##Solution 1      1     2     3
##Solution 2      1     2     3
##Solution 3      1     2     5
##Solution 4      1     2     6
##Solution 5      3     4     6
##Solution 6      3     4     5
##Solution 7      3     4     5
##Solution 8      1     3     6
##Solution 9      2     4     5
##Solution 10     1     3     4
##
##$values
## Solution 1  Solution 2  Solution 3  Solution 4  Solution 5  Solution 6 
##  0.9141995   0.9141995   0.9098502   0.9074543   0.9034868   0.9020271 
## Solution 7  Solution 8  Solution 9 Solution 10 
##  0.9020271   0.8988192   0.8982510   0.8940945 
##
##$bestvalues
##   Card.3 
##0.9141995 
##
##$bestsets
##Var.1 Var.2 Var.3 
##    1     2     3 
}

\keyword{manip}
