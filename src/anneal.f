                
!#######################################################################
       SUBROUTINE anneal(criterio,p,svector,kmin,kmax,valores,vars,
     +   bestval,bestvar,nfora,fora,ndentro,dentro,nsol,
     +   niter,improvement,printfile,iseed,beta,citer,nqsi,qsi,esp)
!#######################################################################

! Routine to perform a Simulated Annealing search for a k-variable subset
! of a set of p original variables, which maximizes one of several criteria.
! (The criteria currently considered are RM, RV and GCD). The user may force the inclusion
! and/or exclusion of certain variables from the k-subset.
! 
! Designed to be called by an R/S function "anneal" (included in the package "subselect",
! available from CRAN), it can be called directly
!  as a Fortran subroutine, provided its arguments are correctly defined (all
! checks and most initializations are in the R/S function).
!
! WARNING: Requires the LAPACK routines SLARUV, DSYEVD and DPOSV
!
! INPUT : 
!
!criterio  - integer variable indicating which criterion of subset quality
!            was requested.
!   p      - integer variable indicating the number of original variables.
! svector  - double precision vector, giving the full covariance (or 
!            correlation) matrix of the p variables (given as a vector
!            for convenience in passing from R to Fortran). 
!kmin,kmax - integer variables, giving the smallest and largest cardinalities
!           of the subsets of variables that are wanted.
!  nfora   - integer variable indicating the number of original variables 
!            that are to be forcefully excluded from the k-subset.
!   fora   - integer vector, indicating the variable numbers associated
!            with the nfora variables that are to be excluded by force.
! ndentro  - integer variable indicating the number of original variables 
!            that are to be forcefully included in the k-subset.
!  dentro  - integer vector, indicating the variable numbers associated
!            with the ndentro variables that are to be included by force.
!  nsol    - integer indicating how many different solutions are to be 
!            produced by the simulated annealing search.
! niter    - integer variable, number of simulated annealing iterations 
!           requested for each solution.
!improvement- logical variable. If true, the best solution produced by 
!            Simulated Annealing will be subjected to a restricted local
!            improvement algorithm. (By default .true.).
!printfile  - logical variable, if .true. the solutions are printed out to a file
!                 (see below for details).
! iseed    - integer vector of length 4, used as a random seed by SLARUV.
!  beta    - double precision variable,(between 0 and 1) indicating the
!           geometric cooling factor for simulated annealing.
! citer    - double precision variable, indicating the initial temperature
!           to be used in the simulated annealing algorithm.
!  nqsi    - integer variable, indicating the dimension of the Principal 
!           Subspace with which the subspace spanned by the k-subset is
!            to be compared. (For GCD criterion only).
!  qsi     - integer vector, giving the ranks of the Principal Components 
!           which span the Principal Subspace described above (see nqsi).
!           (For GCD criterion only).
!  esp     - logical vector. If .true., indicates that the user has specified
!            a set of Principal Components different from the first k, to 
!            compare with the k-variable subset. (For GCD criterion only).
!
! OUTPUT: 
!
! valores  - double precision variable giving the final values produced by S. Annealing 
!                  for each of the nsol solutions, in each cardinality, for the
!                  criterion required. Will pass to the R/S function "anneal"
!                  for conversion to an nsol x (kmax-kmin+1) matrix within R/S, 
!                  with each column associated with a  cardinality and each row with a different 
!                  solution. 
!  vars    - integer vector giving the list of variable numbers  belonging to  each of the nsol
!                final subsets, for each cardinality (padded with zeros if k<kmax). Will pass to the 
!                R/S function "anneal"  for conversion to an nsol x kmax x (kmax-kmin+1) 
!                3-dimensional array.
! bestval - double precision vector with the best value of the criterion
!                 obtained by running the routine, for each cardinality.  If improvement=.true. (the 
!                 default), then the value given is the value resulting from running the restricted local
!                 search algorithm on the best solution from Simulated Annealing (output for R/S).
! bestvar - integer vector with the variable numbers of the best subset obtained, for
!                each cardinality.  If improvement=.true. (the 
!                 default), then the subset given results from running the restricted local
!                 search algorithm on the best solution from Simulated Annealing (output for R/S).
! 
! If printfile = .true. (default is printfile = .false.) output is written, directly from the Fortran code, 
!                 to files, giving the  best values of the criteria, and respective subsets,  for each 
!                 cardinality, as well as indicating possible improvements, if improvement=.true. .
!                  The generations in which new best values are obtained are also indicated, to help
!                  assess whether the genetic algorithm has been running long enough.
!                  The files where these results are written out to (in the working directory) are called:
!                       rmsa.txt  -- RM criterion
!                       rvsa.txt  -- RV criterion
!                       gcdsa.txt -- GCD criterion

! general declarations
       CHARACTER*3 critname(10)
       INTEGER p,iseed(4),criterio,niter
       INTEGER fora(0:nfora),fica(0:p),dentro(0:ndentro),auxw(kmax)
       INTEGER vars(nsol*(kmax-kmin+1)*kmax)
       INTEGER bestvar((kmax-kmin+1)*kmax)
       LOGICAL setk(300),setkmax(300),printfile,improvement
       DOUBLE PRECISION s(300,300),sq(300,300),svector(p*p)
       DOUBLE PRECISION vmax,vactual,vcorrente
       DOUBLE PRECISION valores((kmax-kmin+1)*nsol)
       DOUBLE PRECISION beta,citer,critvalue
       DOUBLE PRECISION bestval(kmax-kmin+1)
! declarations only for the RM criterion
       DOUBLE PRECISION tracos,dobjrm
! declarations only for the RV criterion
       DOUBLE PRECISION tracosq,dobjrv
! declarations only for the GCD criterion
       DOUBLE PRECISION valp(p), vecp(300,300),dobjgcd
       INTEGER nqsi,qsi(p)
       LOGICAL esp

       external randsk1,dobjrm,dobjrv,dobjgcd
       external dcorrigesk,dannealing
       external dmelhoramentogen,inicializar

! initializations

      call inicializar(criterio,p,s,svector,sq,nfora,fora,ndentro,
     +  dentro,fica,tracos,tracosq,valp,vecp)

      if (printfile) then
        critname(1)='RM'
        critname(2)='RV'
        critname(3)='GCD'

        if (criterio.eq.1) then
          open(3,file='rmsa.txt')
        end if      
        if (criterio.eq.2) then
          open(3,file='rvsa.txt')
        end if      
        if (criterio.eq.3) then
          open(3,file='gcdsa.txt')
        end if
        write(3,*) '        Simulated Annealing '
        write(3,*) '         Criterion:  ',critname(criterio)
        write(3,*) ' '
        write(3,*) 'Number of iterations for each solution: ',niter
        write(3,*) ' '
       end if
!*********************************
! The loop which is to be repeated for each cardinality of subsets, 
! between kmin e kmax starts here and ends together with the subroutine.

      do k=kmin,kmax
         if (printfile) then
           write(3,32) ' CARDINALITY   k=',k
 32        format(' ',/,15x,a20,1X,i3,/)  
         end if

! Generates nsol solutions via Simulated Annealing

         if (criterio.eq.3) then
             if (.not.esp) then
               nqsi=k
               do m=1,nqsi
                 qsi(m)=m
               end do
             end if
         end if
         vmax=0

         do ksol=1,nsol
           call randsk1(iseed,p-ndentro,k-ndentro,setk)
           if(ndentro>0) call dcorrigesk(ndentro,dentro,p,setk)
           if (criterio.eq.1) then
               vactual=dobjrm(k,setk,p,s,sq)
           end if
           if (criterio.eq.2) then
               vactual=dobjrv(k,setk,p,s,sq)
           end if
           if (criterio.eq.3) then
               vactual=dobjgcd(nqsi,qsi,valp,vecp,k,setk,p,s,fica)
           end if

           call dannealing(criterio,p,k,s,sq,setk,vactual,ndentro,
     +       dentro,iseed,niter,beta,citer,tracos,tracosq,
     +       nqsi,qsi,valp,vecp,fica)

           if (criterio.eq.1) then
              critvalue = dsqrt(vactual/tracos)
           end if
           if (criterio.eq.2) then
              critvalue = dsqrt(vactual/tracosq)
           end if
           if (criterio.eq.3) then
              critvalue = vactual/dsqrt(dfloat(nqsi*k))
           end if
           
           if (printfile) then
             write(3,*) 'For solution',ksol, '  ',critname(criterio),
     +         ' =',critvalue
           end if

! preparing the output for R/S
            valores((k-kmin)*nsol+ksol)=critvalue
            jjaux=0
           do i=1,p
             if (setk(i)) then
               jjaux=jjaux+1
               vars(nsol*kmax*(k-kmin)+nsol*(jjaux-1)+ksol) = fica(i)
             end if
           end do  

           if(vactual>vmax) then
              vmax=vactual
              do j=1,p
                 setkmax(j)=setk(j)
              end do
           end if
         end do
         naux=0
         do j=1,p
           if(setkmax(j)) then
              naux=naux+1
              auxw(naux)= fica(j)
           end if
         end do
         if (criterio.eq.1) then
            critvalue=dsqrt(vmax/tracos)
         end if
         if (criterio.eq.2) then
            critvalue=dsqrt(vmax/tracosq)
         end if
         if (criterio.eq.3) then
              critvalue = vmax/dsqrt(dfloat(nqsi*k))
         end if
          if (printfile) then 
            write(3,*) ' '
            write(3,24) 'Best ',critname(criterio), ' =',
     +        critvalue,'With subset:',(auxw(j),j=1,naux)    
 24        format(' ',a5,a3,a2,1X,f10.7,/,1X,a12,/,18(i4))
          end if


! If improvement is "true" runs a restricted local improvement algorithm
! on the best of the nsol solutions produced by Simulated Annealing.  

         if (improvement) then
            vcorrente = vmax
            call dmelhoramentogen(criterio,p,setkmax,vmax,ndentro,
     +        dentro,k,s,sq,nqsi,qsi,valp,vecp,fica)
            if (vmax > vcorrente) then 
                   naux=0
                   do j=1,p
                     if(setkmax(j)) then
                       naux=naux+1
                       auxw(naux)= fica(j)
                     end if
                   end do
                   if (criterio.eq.1) then
                       critvalue=dsqrt(vmax/tracos)
                   end if
                   if (criterio.eq.2) then
                       critvalue=dsqrt(vmax/tracosq)
                   end if
                   if (criterio.eq.3) then
                       critvalue = vmax/dsqrt(dfloat(nqsi*k))
                   end if
                   if (printfile) then
                    write(3,25) ' Improved ',critname(criterio),' =',
     +               critvalue, 'Improved subset:',
     +               (auxw(j),j=1,naux)
 25               format(' ',/,a10,a3,a2,1X,f10.7,/,1X,a16,/,18(i4)) 
                   end if
            else
                   if (printfile) then
                        write(3,*) ' '              
                        write(3,*) 'SA solution was not improved'
                   end if
            end if      
         end if

! End of restricted improvement

         bestval(k-kmin+1)=critvalue
         do j=1,k
            bestvar(kmax*(k-kmin)+j)=auxw(j)
         end do
       end do  

! End of the k-loop

      if (printfile) close(3)
      return
      end





