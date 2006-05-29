\name{leaps}

\alias{leaps}
\concept{Leaps and Bounds algorithm}

\title{A Leaps and Bounds Algorithm for finding the best variable subsets} 

\description{An exact Algorithm for optimizing criteria that measure the 
quality of k-dimensional variable subsets as approximations to a given
set of variables, or to a set of its Principal Components.}

\details{
For each cardinality k (with k ranging from \code{kmin} to \code{kmax}),
leaps performs a branch and bound search for the best \code{nsol} variable 
subsets according to a specified criterion. Leaps implements Duarte Silva's 
adaptation (references [2] and [3]) of Furnival and Wilson's Leaps and
Bounds Algorithm (reference [4]) for variable selection in Regression
Analysis. If the search is  not completed within a user defined time
limit, leaps exits with a warning message. 

The user may force variables to be included and/or excluded from the
\emph{k}-subsets.

In order to improve computation times, the bulk of computations are
carried out by C++ routines. Further details about the Algorithm can 
be found in references [2] and [3] and in the comments to the C++ code.
A discussion of the criteria considered can be found in References [1]
and [3].

The function checks for ill-conditioning of the input matrix
(specifically, it checks whether the ratio of the input matrix's
smallest and largest eigenvalues is less than \code{tolval}). For an
ill-conditioned input matrix, execution is aborted. The function
\code{\link{trim.matrix}} may be used to obtain a well-conditioned input
matrix.


In a general descriptive (Principal Components Analysis) setting, the
three criteria Rm, Rv and Gcd can be used to select good k-variable
subsets. Arguments \code{H} and \code{r} are not used in this context.
See reference [1]  and the \code{Examples} for a more detailed
discussion. 

In the setting of a multivariate linear model, \eqn{X = A \Psi  +
U}{X = A B + U},
criteria Ccr12, Tau2, Xi2 and Zeta2 can be used to select subsets
according to their  contribution to an effect characterized by the
violation of  a reference hypothesis, \eqn{C \Psi = 0}{CB = 0} (see
reference [2] for 
further details). In this setting, arguments \code{mat} and \code{H}
should be set respectively to the usual Total (Hypothesis + Error) and
Hypothesis, Sum of Squares and Cross-Products (SSCP) matrices given by: 
\eqn{\code{mat}  = X' (I - P_{\omega}) X }{mat  = X' (I - P_w) X }
and
\eqn{\code{H}  = X' (P_{\Omega} - P_{\omega}) X}{H  = X' (P_O - P_w) X}
where
\eqn{P_{\Omega} =  A(A'A)^- A'}{P_O =  A(A'A)^- A'}
and
\deqn{P_{\omega} =  A(A'A)^- A' - A(A'A)^- C'[C(A'A)^- C']^- C(A'A)^-
A'}{P_w =  A(A'A)^- A' - A(A'A)^- C'[C(A'A)^- C']^- C(A'A)^-
A'}
are projection matrices on the spaces spanned by the columns of A
(space \eqn{\Omega}{O}) and by the linear combinations of these columns that
satisfy the reference hypothesis (space \eqn{\omega}{w}). In these formulae,
\eqn{M'} denotes the transpose of \eqn{M} and
\eqn{M^-} a generalized inverse.   
Argument \code{r} should be set to the expected rank of \code{H},
which equals the difference between the dimensions of \eqn{\Omega}{O} and
\eqn{\omega}{w}. Currently, for reasons of computational efficiency,
criterion Ccr12 is available only when \eqn{\code{r} \leq  3}{\code{r}<=3}. 
Particular cases in this setting include Linear Discriminant Analyis
(LDA), Linear Regression Analysis (LRA), Canonical Correlation
Analysis (CCA) with one set of variables fixed  and several extensions of
these and other classical multivariate methodologies. For LDA, LRA and CCA,
the reference hypothesis merely states that all observation vectors
should have the same expected value and the mat matrix reduces to the
SSCP of the deviations from their global mean. In this case, \code{mat} can
be converted into a covariance matrix if \code{H} is likewise divided by the
same degrees of freedom. 
}

\usage{leaps(mat,kmin=1,kmax=ncol(mat)-1,nsol=1,exclude=NULL,
include=NULL,criterion="default",pcindices="first_k",timelimit=15,H=NULL,
r=0, tolval=10*.Machine$double.eps,tolsym=1000*.Machine$double.eps)} 


\arguments{

  \item{mat}{a covariance/correlation or sums of squares and products
    matrix of the variables from which the k-subset is to be selected. See the \code{Details} section below.}

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
  is to be used in judging the quality of the subsets. Currently, 
  the "Rm", "Rv", "Gcd", "Tau2", "Xi2", "Zeta2" and "Ccr12" criteria are
  supported  (see the \code{Details} section, the \code{References} and
  the links   \code{\link{rm.coef}},  
  \code{\link{rv.coef}}, \code{\link{gcd.coef}}, \code{\link{tau2.coef}},
  \code{\link{xi2.coef}}, \code{\link{zeta2.coef}} and
 \code{\link{ccr12.coef}} for further details).}

  \item{pcindices}{either a vector of ranks of Principal Components that are to be
  used for comparison with the k-variable subsets (for the Gcd
  criterion only, see \code{\link{gcd.coef}}) or the default text
  \code{first_k}. The latter will associate PCs 1 to \emph{k} with each
  cardinality \emph{k} that has been requested by the user.}

  \item{timelimit}{a user specified limit (in seconds) for the maximum time 
   allowed to conduct the search. After this limit is exceeded, leaps exits
   with a waring message stating that it was not possible to find the otpimal
   subsets within the allocated time.}
 
  \item{H}{Effect description matrix. Not used with the Rm, Rv or Gcd
    criteria, hence the NULL default value. See the \code{Details} section below.}
  
  \item{r}{Expected rank of the effects (\code{H}) matrix. Not used with the
    Rm, Rv or Gcd criteria. See the \code{Details} section below.}

 \item{tolval}{the tolerance level for the reciprocal of the 2-norm
condition number of the correlation/covariance or sums of squares matrix, i.e., for the ratio of the smallest to the largest eigenvalue of the input matrix. Matrices with a reciprocal of the condition number smaller than \code{tolval} will abort the search algorithm.}

 \item{tolsym}{the tolerance level for symmetry of the
   covariance/correlation/total matrix and for the effects (\code{H})
   matrix. If corresponding matrix entries differ by more than this
   value, the input matrices will be considered asymmetric and execution
   will be aborted. If corresponding entries are different, but by less
   than this value, the input matrix will be replaced by its symmetric
   part, i.e., input matrix A becomes (A+t(A))/2.}

 }

\value{A list with five items:

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

   \item{call}{The function call which generated the output.}
}

\seealso{\code{\link{rm.coef}}, \code{\link{rv.coef}},
\code{\link{gcd.coef}},  \code{\link{tau2.coef}}, \code{\link{xi2.coef}},
\code{\link{zeta2.coef}}, \code{\link{ccr12.coef}},
\code{\link{anneal}}, \code{\link{genetic}},  
\code{\link{improve}}, \code{\link{trim.matrix}}.}

\references{
[1] Cadima, J. and Jolliffe, I.T. (2001). Variable Selection and the
Interpretation of Principal Subspaces, \emph{Journal of Agricultural,
Biological and Environmental Statistics}, Vol. 6, 62-79.

[2]  Duarte Silva, A.P. (2001) Efficient Variable Screening for
Multivariate Analysis, Journal of Multivariate Analysis Vol. 76,
35-62. 

[3] Duarte Silva, A.P. (2002) Discarding Variables in a Principal 
Component Analysis: Algorithms for All-Subsets Comparisons,
\emph{Computational Statistics}, Vol. 17, 251-271.

[4] Furnival, G.M. and Wilson, R.W. (1974). Regressions by Leaps and
Bounds, \emph{Technometrics}, Vol. 16, 499-511.

}

\examples{

## --------------------------------------------------------------------

## 
## 1) For illustration of use, a small data set.
## Subsets of variables of all cardinalities are sought using the
## RM criterion.
##

data(swiss)
leaps(cor(swiss),nsol=3, criterion="RM")

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
##$call
##leaps(cor(swiss), nsol = 3, criterion="RM")


## --------------------------------------------------------------------

##
## 2) Asking only for 2- and 3- dimensional subsets and excluding 
## variable number 6.
## 

data(swiss)
leaps(cor(swiss),2,3,exclude=6,nsol=3,criterion="gcd")

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
##
##$call
##leaps(cor(swiss), 2, 3, exclude = 6, nsol = 3, criterion = "gcd")

## --------------------------------------------------------------------

##
## 3) Searching for 2- and 3- dimensional subsets that best approximate 
## the spaces generated by the first three Principal Components
##


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
##
##$call
##leaps(cor(swiss), 2, 3, criterion = "gcd", pcindices = 1:3, nsol = 3)


## --------------------------------------------------------------------

##
## 4) An example of subset selection in the context of Multiple Linear
## Regression. Variable 5 in the Cars93 MASS library data set is
## regressed on 13 other variables. A best subset of linear predictors
## is sought, using the "ZETA_2" criterion which, in the case 
## of a Linear Regression, is merely  the standard Coefficient of
## Determination, R^2 (as are the other three criteria for the
## multivariate linear hypothesis, "XI_2", "CCR1_2" and "TAU_2").
##

library(MASS)
data(Cars93)
Y <- Cars93[,5]
X <- Cars93[,c(7:8,12:15,17:22,25)]
Sy <- var(Y)
Sx <- var(X)
Sxy <- var(X,Y)
Carsmat <- Sx
CarsH <- Sxy \%*\% t(Sxy)/Sy
colnames(X)

##  [1] "MPG.city"           "MPG.highway"        "EngineSize"        
##  [4] "Horsepower"         "RPM"                "Rev.per.mile"      
##  [7] "Fuel.tank.capacity" "Passengers"         "Length"            
## [10] "Wheelbase"          "Width"              "Turn.circle"       
## [13] "Weight"            


leaps(Carsmat, kmin=4, kmax=6, H=CarsH, r=1, crit="zeta2")

## $subsets
## , , Card.4
## 
##            Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Solution 1     4     5    10    11     0     0
## 
## , , Card.5
## 
##            Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Solution 1     4     5    10    11    12     0
## 
## , , Card.6
## 
##            Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Solution 1     4     5     9    10    11    12
## 
## 
## $values
##               card.4    card.5   card.6
## Solution 1 0.7143794 0.7241457 0.731015
## 
## $bestvalues
##    Card.4    Card.5    Card.6 
## 0.7143794 0.7241457 0.7310150 
## 
## $bestsets
##        Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Card.4     4     5    10    11     0     0
## Card.5     4     5    10    11    12     0
## Card.6     4     5     9    10    11    12
##

## The selected best variable subsets

colnames(X)[c(4,5,10,11)]

## [1] "Horsepower" "RPM"        "Wheelbase"  "Width"     

colnames(X)[c(4,5,10,11,12)]

## [1] "Horsepower"  "RPM"         "Wheelbase"   "Width"       "Turn.circle"

colnames(X)[c(4,5,9,10,11,12)]

## [1] "Horsepower"  "RPM"         "Length"      "Wheelbase"   "Width"      
## [6] "Turn.circle"


## --------------------------------------------------------------------

## 5) A Linear Discriminant Analysis example with a very small data set. 
## We consider the Iris data and three groups, defined by species (setosa, 
## versicolor and virginica). The goal is to select the 2- and 3-variable
## subsets that are optimal for the linear discrimination (as measured 
## by the "CCR1_2" criterion).
## Both the 'iris' and the 'iris3' R objects are used, for convenience in  
## computing the Total matrix (T=irismat) and the Effects matrix (H=irish).
  
data(iris)
data(iris3) 
irismat <-149*var(iris[,1:4])
weights<-c(49,49,49)
irish <- irismat - matrix(weights\%*\% t(apply(iris3,3,var)),4,4)
leaps(149*var(iris[,c(1,2,3,4)]),kmin=2,kmax=3,H=irish,r=2,crit="ccr12")

## $subsets
## , , Card.2
## 
##            Var.1 Var.2 Var.3
## Solution 1     1     3     0
## 
## , , Card.3
## 
##            Var.1 Var.2 Var.3
## Solution 1     2     3     4
## 
## 
## $values
##               card.2   card.3
## Solution 1 0.9589055 0.967897
## 
## $bestvalues
##    Card.2    Card.3 
## 0.9589055 0.9678971 
## 
## $bestsets
##        Var.1 Var.2 Var.3
## Card.2     1     3     0
## Card.3     2     3     4


## --------------------------------------------------------------------

## 6) An example of subset selection in the context of a Canonical
## Correlation Analysis. Two groups of variables within the Cars93
## MASS library data set are compared. The goal is to select 4- to
## 6-variable subsets of the 13-variable 'X' group that are optimal in
## terms of preserving the canonical correlations, according to the
## "ZETA_2" criterion (Warning: the 3-variable 'Y' group is kept
## intact; subset selection is carried out in the 'X' 
## group only).  The 'tolsym' parameter is used to relax the symmetry
## requirements on the effect matrix H which, for numerical reasons,
## is slightly asymmetric. Since corresponding off-diagonal entries of
## matrix H are different, but by less than tolsym, H is replaced  
## by its symmetric part: (H+t(H))/2.

library(MASS)
data(Cars93)
Y <- Cars93[,4:6]
X <- Cars93[,c(7:8,12:15,17:22,25)]
Sy <- var(Y)
Sx <- var(X)
Sxy <- var(X,Y)
Carsmat <- Sx
CarsH <- Sxy \%*\% solve(Sy,t(Sxy))
colnames(X)

##  [1] "MPG.city"           "MPG.highway"        "EngineSize"        
##  [4] "Horsepower"         "RPM"                "Rev.per.mile"      
##  [7] "Fuel.tank.capacity" "Passengers"         "Length"            
## [10] "Wheelbase"          "Width"              "Turn.circle"       
## [13] "Weight"            


leaps(Carsmat, kmin=4, kmax=6, H=CarsH, r=3, crit="zeta2", tolsym=1e-9)

## $subsets
## , , Card.4
## 
##            Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Solution 1     3     4    10    11     0     0
## 
## , , Card.5
## 
##            Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Solution 1     4     5     9    10    11     0
## 
## , , Card.6
## 
##            Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Solution 1     4     5     9    10    11    12
## 
## 
## $values
##               card.4    card.5    card.6
## Solution 1 0.4827353 0.5018922 0.5168627
## 
## $bestvalues
##    Card.4    Card.5    Card.6 
## 0.4827353 0.5018922 0.5168627 
## 
## $bestsets
##        Var.1 Var.2 Var.3 Var.4 Var.5 Var.6
## Card.4     3     4    10    11     0     0
## Card.5     4     5     9    10    11     0
## Card.6     4     5     9    10    11    12
## 
## Warning message:
## 
##  The effect description matrix (H) supplied was slightly asymmetric: 
##  symmetric entries differed by up to 1.20053300634027e-10.
##  (less than the 'tolsym' parameter).
##  The H matrix has been replaced by its symmetric part.
##  in: validnovcrit(mat, criterion, H, r, p, tolval, tolsym) 

## --------------------------------------------------------------------

}

\keyword{manip}

