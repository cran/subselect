\name{gcd.coef}
\alias{gcd.coef}
\title{COMPUTES YANAI'S GCD FOR THE VARIABLE-SUBSET SELECTION PROBLEM}
\description{
Computes Yanai's Generalized Coefficient of Determination for the
similarity of the subspaces spanned by a subset of
variables and a subset of the full data set's Principal
Components.
}
\usage{gcd.coef(mat, indices, pcindices = seq(1:length(indices)))}
\arguments{
  \item{mat}{the full data set's covariance (or correlation) matrix.}
  \item{indices}{a numerical vector giving the indices of the
variables in the subset.} 
  \item{pcindices}{a numerical vector of indices of Principal
  Components. By default, the first k PCs are chosen, where k is
the cardinality of the subset of variables specified by \code{indices}.}
}
\details{
Computes Yanai's Generalized Coefficient of Determination for the
similarity of the subspaces spanned by a subset of
variables (specified by \code{indices}) and a subset of the
full-data set's Principal Components (specified by \code{pcindices}).
 Input data is expected in the form of a (co)variance or
correlation matrix. If a non-square matrix is given, it is assumed to
be a data matrix, and its (co)variance matrix is used as input. The
number of variables (k) and of PCs (q) does not have to be the same.

Yanai's GCD is defined as:
\deqn{GCD = \frac{\mathrm{tr}(P_v\cdot P_c)}{\sqrt{k\cdot q}}}{GCD = tr(Pv
Pc)/sqrt(k q)}
where \eqn{P_v}{Pv} and \eqn{P_c}{Pc} are the matrices of orthogonal
projections on the subspaces spanned by the k-variable subset and by
the q-Principal Component subset, respectively.

This definition is equivalent to:
\deqn{GCD = \frac{1}{\sqrt{k q}} \sum\limits_{i}(r_m)_i^2}{GCD = sum_i
(r_i^2) / sqrt(k q)} 
where \eqn{(r_m)_i}{r_i} stands for the multiple correlation between the
\code{i}-th Principal Component and the k-variable subset, and the sum
is carried out over the q PCs (i=1,...,q) selected.

These definitions are also equivalent to the expression used in the
code, which only requires the covariance (or correlation) matrix of
the data under consideration.
}
\value{
The value of the GCD coefficient.
}
\references{
 Cadima, J. and Jolliffe, I.T. (2001), "Variable Selection and the
Interpretation of Principal Subspaces", \emph{Journal of Agricultural,
Biological and Environmental Statistics}, Vol. 6, 62-79.

 Ramsay, J.O., ten Berge, J. and Styan, G.P.H. (1984), "Matrix
Correlation", \emph{Psychometrika}, 49, 403-423. 
}
\examples{
data(iris3) 
x<-iris3[,,1]
gcd.coef(cor(x),c(1,3))
## [1] 0.7666286
gcd.coef(cor(x),c(1,3),pcindices=c(1,3))
## [1] 0.584452
gcd.coef(cor(x),c(1,3),pcindices=1)
## [1] 0.6035127
}
\keyword{manip}
