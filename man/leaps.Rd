\name{leaps}

\alias{leaps}

\title{A Leaps and Bounds Algorithm for finding the best variable subsets} 

\description{An exact Algorithm for optimizing criteria that measure the 
quality of k-dimensional variable subsets as approximations to a given
set of variables, or to a set of its Principal Components.}

\details{
For each cardinality k (with k ranging from \code{kmin} to \code{kmax}),
leaps performs a branch and bound search for the best \code{nsol} variable 
subsets according to a specified criterion. Leaps implements Duarte Silva's 
adaptation (Reference 2) of Furnival and Wilson's Leaps and Bounds Algorithm  
(Reference 3) for variable selection in Regression Analysis. If the search is 
not completed within an user defined time limit,leaps exits with a warning 
message. 

The user may force variables to be included and/or excluded from the
\emph{k}-subsets.

In order to improve computation times, the bulk of computations are
carried out by C++ routines. Further details about the Algorithm can 
be found in References 2 and 3.
A discussion of the criteria considered can be found in Reference 1.}

\usage{leaps(mat,kmin=1,kmax=ncol(mat)-1,nsol=1,exclude=NULL,include=NULL,
criterion="RM",pcindices=1:kmin,timelimit=15)}


\arguments{
  \item{mat}{a covariance or correlation matrix of the variables from
  which the k-subset is to be selected.}

  \item{kmin}{the cardinality of the smallest subset that is wanted.}

  \item{kmax}{the cardinality of the largest subset that is wanted.}

  \item{nsol}{the number of different subsets of each cardinality that
   are requested .}

  \item{exclude}{a vector of variables (referenced by their row/column
  numbers in matrix \code{mat}) that are to be forcibly excluded from
  the subsets.} 

  \item{include}{a vector of variables (referenced by their row/column
  numbers in matrix \code{mat}) that are to be forcibly included in
  the subsets.} 

  \item{criterion}{Character variable, which indicates which criterion
  is to be used in judging the quality of the subsets. Currently, only
  the RM, RV and GCD criteria are supported, and referenced as "RM",
  "RV" or "GCD" (see References, \code{\link{rm.coef}}, 
  \code{\link{rv.coef}} and \code{\link{gcd.coef}} for further details).}

  \item{pcindices}{a vector of ranks of Principal Components that are to be
  used for comparison with the variable subsets (for the GCD
  criterion only, see \code{\link{gcd.coef}}).}


  \item{timelimit}{a user specified limit (in seconds) for the maximum time 
   allowed to conduct the search. After this limit is exceeded, leaps exits
   with a waring message stating that it was not possible to find the otpimal
   subsets within the allocated time.}
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
   of the best \code{nsol} (rows) subsets according to the chosen criterion.} 

   \item{bestvalues}{A length(\code{kmin}:\code{kmax}) vector giving
   the overall best values of the criterion for each cardinality.} 

   \item{bestsets}{A length(\code{kmin}:\code{kmax}) x \code{kmax}
   matrix, giving, for each cardinality (rows), the variables
   (referenced by their row/column numbers in matrix \code{mat}) in
   the best k-subset.}
}

\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}}, \code{\link{anneal}}, \code{\link{genetic}}, 
\code{\link{improve}}.}

\references{
1) Cadima, J. and Jolliffe, I.T. (2001). Variable Selection and the
Interpretation of Principal Subspaces, \emph{Journal of Agricultural,
Biological and Environmental Statistics}, Vol. 6, 62-79.

2) Duarte Silva, A.P. (2002) Discarding Variables in a Principal 
Component Analysis: Algorithms for All-Subsets Comparisons,
\emph{Computational Statistics}, Vol. 17, 251-271.

3) Furnival, G.M. and Wilson, R.W. (1974). Regressions by Leaps and
Bounds, \emph{Technometrics}, Vol. 16, 499-511.
}

\examples{

# For illustration of use, a small data set with very few iterations
# of the algorithm. 

data(swiss)
leaps(cor(swiss),nsol=3,criterion="RM")

##$subsets
##, , Card.1
##
##           Var.1 Var.2 Var.3 Var.4 Var.5
##Solution 1     3     0     0     0     0
##Solution 2     1     0     0     0     0
##Solution 3     4     0     0     0     0
##
##, , Card.2
##
##           Var.1 Var.2 Var.3 Var.4 Var.5
##Solution 1     3     6     0     0     0
##Solution 2     4     5     0     0     0
##Solution 3     1     2     0     0     0
##
##, , Card.3
##
##           Var.1 Var.2 Var.3 Var.4 Var.5
##Solution 1     4     5     6     0     0
##Solution 2     1     2     5     0     0
##Solution 3     3     4     6     0     0
##
##, , Card.4
##
##           Var.1 Var.2 Var.3 Var.4 Var.5
##Solution 1     2     4     5     6     0
##Solution 2     1     2     5     6     0
##Solution 3     1     4     5     6     0
##
##, , Card.5
##
##           Var.1 Var.2 Var.3 Var.4 Var.5
##Solution 1     1     2     3     5     6
##Solution 2     1     2     4     5     6
##Solution 3     2     3     4     5     6
##
##
##$values
##              card.1    card.2    card.3    card.4    card.5
##Solution 1 0.6729689 0.8016409 0.9043760 0.9510757 0.9804629
##Solution 2 0.6286185 0.7982296 0.8791856 0.9506434 0.9776338
##Solution 3 0.6286130 0.7945390 0.8777509 0.9395708 0.9752551
##
##$bestvalues
##   Card.1    Card.2    Card.3    Card.4    Card.5 
##0.6729689 0.8016409 0.9043760 0.9510757 0.9804629 
##
##$bestsets
##       Var.1 Var.2 Var.3 Var.4 Var.5
##Card.1     3     0     0     0     0
##Card.2     3     6     0     0     0
##Card.3     4     5     6     0     0
##Card.4     2     4     5     6     0
##Card.5     1     2     3     5     6
##

#
#
# Asking only for 2- and 3- dimensional subsets and excluding 
# variable number 6.
# 

data(swiss)
leaps(cor(swiss),2,3,exclude=6,nsol=3,criterion="RM")

##$subsets
##, , Card.2
##
##           Var.1 Var.2 Var.3
##Solution 1     4     5     0
##Solution 2     1     2     0
##Solution 3     1     3     0
##
##, , Card.3
##
##           Var.1 Var.2 Var.3
##Solution 1     1     2     5
##Solution 2     1     4     5
##Solution 3     2     4     5
##
##
##$values
##              card.2    card.3
##Solution 1 0.7982296 0.8791856
##Solution 2 0.7945390 0.8686515
##Solution 3 0.7755232 0.8628693
##
##$bestvalues
##   Card.2    Card.3 
##0.7982296 0.8791856 
##
##$bestsets
##       Var.1 Var.2 Var.3
##Card.2     4     5     0
##Card.3     1     2     5

#
# Searching for 2- and 3- dimensional subsets that best approximate the 
# spaces generated by the first three Principal Components
#


data(swiss)
leaps(cor(swiss),2,3,criterion="gcd",pcindices=1:3,nsol=3)

##$subsets
##, , Card.2
##
##           Var.1 Var.2 Var.3
##Solution 1     4     5     0
##Solution 2     5     6     0
##Solution 3     4     6     0
##
##, , Card.3
##
##           Var.1 Var.2 Var.3
##Solution 1     4     5     6
##Solution 2     3     5     6
##Solution 3     2     5     6
##
##
##$values
##              card.2    card.3
##Solution 1 0.7831827 0.9253684
##Solution 2 0.7475630 0.8459302
##Solution 3 0.7383665 0.8243032
##
##$bestvalues
##   Card.2    Card.3 
##0.7831827 0.9253684 
##
##$bestsets
##       Var.1 Var.2 Var.3
##Card.2     4     5     0
##Card.3     4     5     6

}

\keyword{manip}
