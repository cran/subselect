\name{anneal}

\alias{anneal}

\title{Simulated Annealing Search for an optimal k-variable subset} 

\description{Given a set of variables,  a Simulated Annealing
algorithm seeks a k-variable subset which is
optimal, as a surrogate for the whole set, with respect to a given criterion. 
}

\details{An initial k-variable subset (for k ranging from \code{kmin} to \code{kmax}) 
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
algorithm is the number of iterations (\code{niter}) which is controlled by the
user. Also controlled by the user are the initial temperature (\code{temp}) the
rate of geometric cooling of the temperature (\code{cooling}) and the 
frequency with which the temperature is cooled, as measured by 
\code{coolfreq}, the number of iterations after which the temperature is
multiplied by 1-\code{cooling}.

Optionally, the best k-variable subset produced by Simulated Annealing may be
passed as input to a restricted local search algorithm, for possible further
improvement.

The user may force variables to be included and/or excluded from the
k-subsets, and may specify initial solutions.

For each cardinality k, the total number of calls to the procedure
which computes the criterion 
values is \code{nsol} x (\code{niter} + 1). These calls are the
dominant computational effort in each iteration of the algorithm.  

In order to improve computation times, the bulk of computations is
carried out by a Fortran routine. Further details about the Simulated
Annealing algorithm can 
be found in Reference 1 and in the comments to the Fortran code (in
the \code{src} subdirectory for this package). For p > 300, it is 
necessary to change the declarative statements in the Fortran code. 
}

\usage{anneal( mat, kmin, kmax = kmin, nsol = 1, niter = 1000, exclude
= NULL, include = NULL, improvement = TRUE, setseed = FALSE,
cooling = 0.05,  temp = 1, coolfreq = 1, criterion = "RM",
pcindices = 1:kmax), initialsol=c(0) )} 

\arguments{

  \item{mat}{a covariance or correlation matrix of the variables from
  which the k-subset is to be selected.}

  \item{kmin}{the cardinality of the smallest subset that is wanted.}

  \item{kmax}{the cardinality of the largest subset that is wanted.}

  \item{nsol}{the number of initial/final subsets (runs of the algorithm).}

  \item{niter}{the number of iterations of the algorithm for each
  initial subset.}  

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

  \item{cooling}{variable in the ]0,1[ interval indicating the rate of
  geometric cooling for the Simulated Annealing algorithm.}

  \item{temp}{positive variable indicating the initial
  temperature for the Simulated Annealing algorithm.}

  \item{coolfreq}{positive integer indicating the number of iterations
  of the algorithm between coolings of the temperature. By default,
  the temperature is cooled at every iteration.}

  \item{criterion}{Character variable, which indicates which criterion
  is to be used in judging the quality of the subsets. Currently, only
  the RM, RV and GCD criteria are supported, and referenced as "RM",
  "RV" or "GCD" (see References, \code{\link{rm.coef}}, 
  \code{\link{rv.coef}} and \code{\link{gcd.coef}} for further
  details).}

  \item{pcindices}{a vector of ranks of Principal Components that are to be
  used for comparison with the k-variable subsets (for the GCD
  criterion only, see \code{\link{gcd.coef}}).}

  \item{initialsol}{vector, matrix or 3-d array of initial solutions
  for the simulated annealing search. If a \emph{single cardinality} is
  required, \code{initialsol} may be a vector of length \emph{k}, in
  which case it is used as the initial solution for all \code{nsol}
  final solutions that are requested; a 1 x \emph{k} matrix (as
  produced by the \code{$bestsets} output value of the algorithm functions
  \code{anneal}, \code{\link{genetic}}, or \code{\link{improve}}, or
  a 1 x \emph{k} x 1 array (as produced by the
  \code{$subsets} output value), in
  which case it will be treated as the above k-vector; or an
  \code{nsol} x \code{k} matrix, or  \code{nsol} x \code{k} x 1 3-d
  array, in which case each row (dimension 1) will be used 
  as the initial solution for each of the \code{nsol} final solutions
  requested. If \emph{more than one cardinality} is requested,
  \code{initialsol} can be a 
  \code{length(kmin:kmax)} x \code{kmax} matrix (as produced by the
  \code{$bestsets} option of the algorithm functions), in which case
  each row will be replicated to produced the initial solution for all
  \code{nsol} final solutions requested in each cardinality, or a
  \code{nsol} x \code{kmax} x \code{length(kmin:kmax)} 3-d array (as
  produced by the 
  \code{$subsets} output option), in which case each row (dimension 1)
  is interpreted as a different initial solution.

  If the \code{exclude} and/or \code{include} options are used,
  \code{initialsol} must also respect those requirements. }
}

\value{A list with four items:

   \item{subsets}{An \code{nsol} x \code{kmax} x
   length(\code{kmin}:\code{kmax}) 3-dimensional 
   array, giving for each cardinality (dimension 3) and each solution
   (dimension 1) the list of variables (referenced by their row/column
  numbers in matrix \code{mat}) in the subset (dimension 2). (For
  cardinalities  smaller than \code{kmax}, the extra final positions are set
  to zero).} 

   \item{values}{An \code{nsol} x length(\code{kmin}:\code{kmax})
   matrix, giving for each cardinality (columns), the criterion values
   of the \code{nsol} (rows) subsets obtained.}

   \item{bestvalues}{A length(\code{kmin}:\code{kmax}) vector giving
   the best values of the criterion obtained for each cardinality. If
   \code{improvement} is TRUE, these values result from the final
   restricted local search algorithm (and may therefore exceed the
   largest value for that cardinality in \code{values}).}

   \item{bestsets}{A length(\code{kmin}:\code{kmax}) x \code{kmax}
   matrix, giving, for each cardinality (rows), the variables
   (referenced by their row/column numbers in matrix \code{mat}) in
   the best k-subset that was found.}
}

\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}}, \code{\link{genetic}}, \code{\link{improve}}.}

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
anneal(cor(swiss),2,3,nsol=4,niter=10,criterion="RM")

##$subsets
##, , Card.2
##
##           Var.1 Var.2 Var.3
##Solution 1     3     6     0
##Solution 2     4     5     0
##Solution 3     1     2     0
##Solution 4     3     6     0
##
##, , Card.3
##
##           Var.1 Var.2 Var.3
##Solution 1     4     5     6
##Solution 2     3     5     6
##Solution 3     3     4     6
##Solution 4     4     5     6
##
##
##$values
##              card.2    card.3
##Solution 1 0.8016409 0.9043760
##Solution 2 0.7982296 0.8769672
##Solution 3 0.7945390 0.8777509
##Solution 4 0.8016409 0.9043760
##
##$bestvalues
##   Card.2    Card.3 
##0.8016409 0.9043760 
##
##$bestsets
##       Var.1 Var.2 Var.3
##Card.2     3     6     0
##Card.3     4     5     6


#
#
# Excluding variable number 6 from the subsets.
# 

data(swiss)
anneal(cor(swiss),2,3,nsol=4,niter=10,criterion="RM",exclude=c(6))

##$subsets
##, , Card.2
##
##           Var.1 Var.2 Var.3
##Solution 1     4     5     0
##Solution 2     4     5     0
##Solution 3     4     5     0
##Solution 4     4     5     0
##
##, , Card.3
##
##           Var.1 Var.2 Var.3
##Solution 1     1     2     5
##Solution 2     1     2     5
##Solution 3     1     2     5
##Solution 4     1     4     5
##
##
##$values
##              card.2    card.3
##Solution 1 0.7982296 0.8791856
##Solution 2 0.7982296 0.8791856
##Solution 3 0.7982296 0.8791856
##Solution 4 0.7982296 0.8686515
##
##$bestvalues
##   Card.2    Card.3 
##0.7982296 0.8791856 
##
##$bestsets
##       Var.1 Var.2 Var.3
##Card.2     4     5     0
##Card.3     1     2     5


# specifying initial solutions: using the subsets produced by
# simulated annealing for one criterion (RM, by default) as initial
# solutions for the simulated annealing search with a different criterion.

data(swiss)
rmresults<-anneal(cor(swiss),2,3,nsol=4,niter=10, setseed=TRUE)
anneal(cor(swiss),2,3,nsol=4,niter=10,criterion="gcd", initialsol=rmresults$subsets)

##$subsets
##, , Card.2
##
##           Var.1 Var.2 Var.3
##Solution 1     3     6     0
##Solution 2     3     6     0
##Solution 3     3     6     0
##Solution 4     3     6     0
##
##, , Card.3
##
##           Var.1 Var.2 Var.3
##Solution 1     4     5     6
##Solution 2     4     5     6
##Solution 3     3     4     6
##Solution 4     4     5     6
##
##
##$values
##              card.2   card.3
##Solution 1 0.8487026 0.925372
##Solution 2 0.8487026 0.925372
##Solution 3 0.8487026 0.798864
##Solution 4 0.8487026 0.925372
##
##$bestvalues
##   Card.2    Card.3 
##0.8487026 0.9253720 
##
##$bestsets
##       Var.1 Var.2 Var.3
##Card.2     3     6     0
##Card.3     4     5     6
}
\keyword{manip}
