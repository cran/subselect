                
!#######################################################################
       SUBROUTINE anneal(criterio,p,svector,kmin,kmax,valores,vars,
     + bestval,bestvar,nfora,fora,ndentro,dentro,nsol,niter,
     + improvement,beta,temp,coolfreq,nqsi,qsi,esp,silog,solinit)
!#######################################################################

! Routine to perform a Simulated Annealing search for a k-variable subset
! of a set of p original variables, which maximizes one of several criteria.
! (The criteria currently considered are RM, RV and GCD). The user may force
!  the inclusion and/or exclusion of certain variables from the k-subset.
! 
! Designed to be called by an R/S function "anneal" (included in the package 
! "subselect", available from CRAN). 
!
! WARNING: Requires the LAPACK routines DSYEVD and DPOSV
! WARNING: Uses R's default Random Number Generator.
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
!  beta    - double precision variable,(between 0 and 1) indicating the
!           geometric cooling factor for simulated annealing.
! temp    - double precision variable, indicating the initial temperature
!           to be used in the simulated annealing algorithm.
! coolfreq - integer indicating the cooling frequency, i.e., the number 
!           of iterations that are to be carried out for each fixed
!           temperature 
!  nqsi    - integer variable, indicating the dimension of the Principal 
!           Subspace with which the subspace spanned by the k-subset is
!            to be compared. (For GCD criterion only).
!  qsi     - integer vector, giving the ranks of the Principal Components 
!           which span the Principal Subspace described above (see nqsi).
!           (For GCD criterion only).
!  esp     - logical vector. If .true., indicates that the user has specified
!            a set of Principal Components different from the first k, to 
!            compare with the k-variable subset. (For GCD criterion only).
! silog    - logical variable indicating whether or not the user has specified
!             initial solutions.
! solinit  - integer variable, giving the user-specified initial solutions
!             (if any)

!
! OUTPUT: 
!
! valores  - double precision variable giving the final values produced by S.  
!            Annealing for each of the nsol solutions, in each cardinality, 
!            for the criterion required. Will pass to the R/S function "anneal"
!            for conversion to an nsol x (kmax-kmin+1) matrix within R/S, 
!            with each column associated with a  cardinality and each row with
!            a different solution. 
!  vars    - integer vector giving the list of variable numbers  belonging 
!            to  each of the nsol final subsets, for each cardinality (padded
!            with zeros if k<kmax). Will pass to the R/S function "anneal" for
!             conversion to an nsol x kmax x (kmax-kmin+1) 3-dimensional array.
! bestval - double precision vector with the best value of the criterion
!           obtained by running the routine, for each cardinality.  If 
!           improvement=.true. (the default), then the value given is the 
!           value resulting from running the restricted local search algorithm
!           on the best solution from Simulated Annealing (output for R/S).
! bestvar - integer vector with the variable numbers of the best subset 
!           obtained, for each cardinality.  If improvement=.true. (the 
!           default), then the subset given results from running the 
!           restricted local search algorithm on the best solution from 
!           Simulated Annealing (output for R/S).
!

! general declarations
       INTEGER p,criterio,niter,poriginal
       INTEGER fora(0:nfora),fica(0:p),dentro(0:ndentro),auxw(kmax)
       INTEGER vars(nsol*(kmax-kmin+1)*kmax)
       INTEGER bestvar((kmax-kmin+1)*kmax)
       LOGICAL setk(300),setkmax(300),improvement
       DOUBLE PRECISION s(300,300),sq(300,300),svector(p*p)
       DOUBLE PRECISION vmax,vactual,vcorrente
       DOUBLE PRECISION valores((kmax-kmin+1)*nsol)
       DOUBLE PRECISION beta,temp,critvalue
       integer coolfreq 
       DOUBLE PRECISION bestval(kmax-kmin+1)
       LOGICAL silog
       INTEGER solinit(kmax*nsol*(kmax-kmin+1))
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
     +  dentro,fica,tracos,tracosq,valp,vecp,poriginal)


!*********************************
! The loop which is to be repeated for each cardinality of subsets, 
! between kmin e kmax starts here and ends together with the subroutine.

      do k=kmin,kmax

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

! determining the initial solution.

           if (silog) then
             do m=1,poriginal
               setk(m)=.false.
             end do           
             do m=1,k
               setk(solinit(nsol*kmax*(k-kmin)+nsol*(m-1)+ksol))=.true.
             end do
             if (nfora.ge.1) then
               do m=1,p
                 setk(m) = setk(fica(m))
               end do
             end if  
           else
               call randsk1(p-ndentro,k-ndentro,setk)           
               if(ndentro>0) call dcorrigesk(ndentro,dentro,p,setk)
           endif



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
     +       dentro,niter,beta,temp,coolfreq,tracos,tracosq,
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

            end if      
         end if

! End of restricted improvement

         bestval(k-kmin+1)=critvalue
         do j=1,k
            bestvar(kmax*(k-kmin)+j)=auxw(j)
         end do
       end do  

! End of the k-loop

      return
      end





