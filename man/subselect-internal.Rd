\name{subselect-internal}
\title{Internal subselect Objects}

\alias{validannimp}
\alias{validation}
\alias{validgenetic}
\alias{validmat}
\alias{validnovcrit}

\description{
  Internal subselect functions that validate input to the search functions.
}
\usage{

validannimp(kmin, kmax, nsol, exclude, nexclude, include, ninclude,
initialsol,implog)

validation(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval,tolsym)

validgenetic(kmin, kmax, popsize, mutprob, exclude, nexclude, include, ninclude, initialpop)

validmat(mat,p,tolval,tolsym)

validnovcrit(mat,criterion,H,r,p,tolval,tolsym)
}
\details{
  These functions are not to be called by the user.
}
\keyword{internal}
