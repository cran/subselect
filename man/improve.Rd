\name{improve}

\alias{improve}

\title{Restricted Local Improvement search for an optimal k-variable subset} 

\description{Given a set of variables,  a Restricted Local Improvement
algorithm seeks a k-variable subset which is
optimal, as a surrogate for the whole set, with respect to a given criterion. 
}

\details{An initial k-variable subset (for k ranging from \code{kmin}
to \code{kmax})  of a full set of p (p not exceeding 300) variables is
randomly selected and the variables not belonging to this subset are
placed in a queue. The possibility of replacing a variable in the
current k-subset with a variable from the queue is then explored.  
More precisely, a variable is selected, removed from the queue, and the 
k values of the criterion which would result from swapping this selected
variable with each variable in the current subset are computed. If the
best of these values improves the current criterion value, the current
subset is updated accordingly. In this case, the variable which leaves
the subset is added to the queue, but only if it has not previously
been in the queue (i.e., no variable can enter the queue twice). The
algorithm proceeds until the queue is emptied.  

The user may force variables to be included and/or excluded from the
k-subsets, and may specify initial solutions.

For each cardinality k, the total number of calls to the procedure
which computes the criterion 
values is O(\code{nsol} x k x p). These calls are the dominant computational
effort in each iteration of the algorithm.  

In order to improve computation times, the bulk of computations are
carried out in a Fortran routine. Further details about the algorithm can
be found in Reference 1 and in the comments to the Fortran code (in
the \code{src} subdirectory for this package). For p > 300, it is
necessary to change the declarative statements in the Fortran code. 
}

\usage{improve( mat, kmin, kmax = kmin, nsol = 1, exclude = NULL,
include = NULL, setseed = FALSE, criterion = "RM", pcindices
= 1:kmax, initialsol=NULL)
}

\arguments{

  \item{mat}{a covariance or correlation matrix of the variables from
  which the k-subset is to be selected.}

  \item{kmin}{the cardinality of the smallest subset that is wanted.}

  \item{kmax}{the cardinality of the largest subset that is wanted.}

  \item{nsol}{the number of different subsets (runs of the algorithm) wanted.}

  \item{exclude}{a vector of variables (referenced by their row/column
  numbers in matrix \code{mat}) that are to be forcibly excluded from
  the subsets.} 

  \item{include}{a vector of variables (referenced by their row/column
  numbers in matrix \code{mat}) that are to be forcibly included from
  the subsets.} 

  \item{setseed}{logical variable indicating whether to fix an initial 
  seed for the random number generator, which will be re-used in future
  calls to this function whenever setseed is again set to TRUE.}

  \item{criterion}{Character variable, which indicates which criterion
  is to be used in judging the quality of the subsets. Currently, only
  the RM, RV and GCD criteria are supported, and referenced as "RM",
  "RV" or "GCD" (see References, \code{\link{rm.coef}}, 
  \code{\link{rv.coef}} and \code{\link{gcd.coef}} for further details).}

  \item{pcindices}{a vector of ranks of Principal Components that are to be
  used for comparison with the k-variable subsets (for the GCD
  criterion only, see \code{\link{gcd.coef}}).}

  \item{initialsol}{vector, matrix or 3-d array of initial solutions
  for the restricted local improvement search. If a \emph{single
  cardinality} is 
  required, \code{initialsol} may be a vector of length
  \emph{k}(accepted even if \code{nsol} > 1, in
  which case it is used as the initial solution for all \code{nsol}
  final solutions that are requested with a warning that the same
  initial solution necessarily produces the same final solution); 
  a 1 x \emph{k} matrix (as
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
  \code{$bestsets} option of the algorithm functions) (even if
  \code{nsol} > 1, in which case
  each row will be replicated to produced the initial solution for all
  \code{nsol} final solutions requested in each cardinality, with a
  warning that a single initial solution necessarily produces
  identical final solutions), or a
  \code{nsol} x \code{kmax} x \code{length(kmin:kmax)} 3-d array (as
  produced by the 
  \code{$subsets} output option), in which case each row (dimension 1)
  is interpreted as a different initial solution.

  If the \code{exclude} and/or \code{include} options are used,
  \code{initialsol} must also respect those requirements. }
}

\value{A list with four items:

   \item{subsets}{An \code{nsol} x \code{kmax} x
   length(\code{kmin}:\code{kmax}) 3-dimensional array, giving for
   each cardinality (dimension 3) and each solution (dimension 1) the
   list of variables (referenced by their row/column numbers in matrix
   \code{mat}) in the subset (dimension 2). (For cardinalities
   smaller than \code{kmax}, the extra final positions are set to zero).}

   \item{values}{An \code{nsol} x length(\code{kmin}:\code{kmax})
   matrix, giving for each cardinality (columns), the criterion values
   of the \code{nsol} (rows) solutions obtained.}

   \item{bestvalues}{A length(\code{kmin}:\code{kmax}) vector giving
   the best values of the criterion obtained for each cardinality.}

   \item{bestsets}{A length(\code{kmin}:\code{kmax}) x \code{kmax}
   matrix, giving, for each cardinality (rows), the variables
   (referenced by their row/column numbers in matrix \code{mat}) in the
   best k-subset that was found.}
}

\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}}, \code{\link{genetic}}, \code{\link{anneal}}.}

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
#

data(swiss)
improve(cor(swiss),2,3,nsol=4,criterion="GCD")
## $subsets
## , , Card.2
##
##            Var.1 Var.2 Var.3
## Solution 1     3     6     0
## Solution 2     3     6     0
## Solution 3     3     6     0
## Solution 4     3     6     0
##
## , , Card.3
##
##            Var.1 Var.2 Var.3
## Solution 1     4     5     6
## Solution 2     4     5     6
## Solution 3     4     5     6
## Solution 4     4     5     6
##
##
## $values
##               card.2   card.3
## Solution 1 0.8487026 0.925372
## Solution 2 0.8487026 0.925372
## Solution 3 0.8487026 0.925372
## Solution 4 0.8487026 0.925372
##
## $bestvalues
##    Card.2    Card.3 
## 0.8487026 0.9253720 
##
## $bestsets
##        Var.1 Var.2 Var.3
## Card.2     3     6     0
## Card.3     4     5     6

#
#
# Forcing the inclusion of variable 1 in the subset
#

 improve(cor(swiss),2,3,nsol=4,criterion="GCD",include=c(1))

## $subsets
## , , Card.2
##
##            Var.1 Var.2 Var.3
## Solution 1     1     6     0
## Solution 2     1     6     0
## Solution 3     1     6     0
## Solution 4     1     6     0
##
## , , Card.3
##
##            Var.1 Var.2 Var.3
## Solution 1     1     5     6
## Solution 2     1     5     6
## Solution 3     1     5     6
## Solution 4     1     5     6
##
##
## $values
##               card.2    card.3
## Solution 1 0.7284477 0.8048528
## Solution 2 0.7284477 0.8048528
## Solution 3 0.7284477 0.8048528
## Solution 4 0.7284477 0.8048528
##
## $bestvalues
##    Card.2    Card.3 
## 0.7284477 0.8048528 
##
## $bestsets
##        Var.1 Var.2 Var.3
## Card.2     1     6     0
## Card.3     1     5     6
}
\keyword{manip}
