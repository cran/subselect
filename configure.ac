AC_PREREQ(2.50)
AC_INIT([src/Makevars.in])
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME], 1)
fi
CC=`grep "^CC" ${R_HOME}/etc/Makeconf | sed "s/^CC *= *//"`
AC_PROG_CC
F77=`grep "^FC" ${R_HOME}/etc/Makeconf | sed "s/^FC *= *//"`
if test -z "${F77}"; then
  AC_MSG_ERROR([R was configured without a FORTRAN compiler], 1)
fi
AC_PROG_F77
FLIBS=`grep "^FLIBS" ${R_HOME}/etc/Makeconf | sed "s/^FLIBS *= *//"`
ACX_BLAS([], AC_MSG_ERROR([could not find your BLAS library], 1))
ACX_LAPACK([], AC_MSG_ERROR([could not find your LAPACK library], 1))
AC_OUTPUT([src/Makevars])
