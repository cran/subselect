                
!#######################################################################
       SUBROUTINE genetic(criterio,p,svector,kmin,kmax,valores,vars,
     +   bestval,bestvar,nfora,fora,ndentro,dentro,npopinic,nger,nclone,
     +   melhorar,improvement,printfile,iseed,nqsi,qsi,esp,kabort)
!#######################################################################

! Routine to perform a Genetic Algorithm search for a k-variable subset
! of a set of p original variables, which maximizes one of several criteria.
! (The criteria currently considered are RM, RV and GCD). The user may force the inclusion
! and/or exclusion of certain variables from the k-subset.
! 
! Designed to be called by an R/S function "genetic" (included in the package "subselect",
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
! npopinic  - integer indicating the size of the initial population of k-subsets.
! nger    - integer variable, number of generations for the genetic algorithm. 
! nclone  - integer variable indicating the maximum number of clones (identical
!                 individuals) that will be tolerated. Clones in excess of nclone are excluded.
! melhorar  - logical variable indicating whether each new child is to be mutated, via
!                a restricted local improvement algorithm, before undergoing survival process.
!improvement- logical variable. If true, the best solution produced by 
!            the Genetic Algorithm will be subjected to a restricted local
!            improvement algorithm. (By default .true.).
!printfile  - logical variable, if .true. the solutions for all iterations 
!           are printed out to a file (see below for details).
! iseed    - integer vector of length 4, used as a random seed by SLARUV.
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
! valores  - double precision variable giving the final values produced by the Genetic Algorithm for 
!                  each of the npopinic members of the population, in each cardinality, for the
!                  criterion required. Will pass to the R/S function "genetic"
!                  for conversion to an npopinic x (kmax-kmin+1) matrix within R/S, 
!                  with each column associated with a  cardinality and each row with a different 
!                  solution. 
!  vars    - integer vector giving the list of variable numbers  belonging to  each of the npopinic
!                final subsets, for each cardinality (padded with zeros if k<kmax). Will pass to the 
!                R/S function "genetic"  for conversion to an npopinic x kmax x (kmax-kmin+1) 
!                3-dimensional array.
! bestval - double precision vector with the best value of the criterion
!                 obtained by running the routine, for each cardinality.  If improvement=.true. (the 
!                 default), then the value given is the value resulting from running the restricted local
!                 search algorithm on the best solution from the Genetic Algorithm (output for R/S).
! bestvar - integer vector with the variable numbers of the best subset obtained, for
!                each cardinality.  If improvement=.true. (the 
!                 default), then the subset given results from running the restricted local
!                 search algorithm on the best solution from the Genetic Algorithm (output for R/S).
! kabort - integer variable indicating the cardinality in which an execution of the
!                algorithm may have been aborted for lack of genetic variability (ouput for R/S).
! 
! If printfile = .true. (default is printfile = .false.) output is written, directly from the Fortran code, 
!                 to files, giving the  best values of the criteria, and respective subsets,  for each 
!                 cardinality, as well as indicating possible improvements, if improvement=.true. .
!                  The generations in which new best values are obtained are also indicated, to help
!                  assess whether the genetic algorithm has been running long enough.
!                  The files where these results are written out to (in the working directory) are called:
!                        rmgen.txt  -- RM criterion
!                        rvgen.txt  -- RV criterion
!                        gcdgen.txt -- GCD criterion


! general declarations
       CHARACTER*3 critname(10)
       INTEGER p,kmin,kmax,iseed(4),criterio,npopinic,nger,nclone
       INTEGER fora(0:nfora),fica(0:p),dentro(0:ndentro),auxw(kmax)
       INTEGER talvez,painaomae(kmax),maenaopai(kmax)
       INTEGER ordemger(0:npopinic)
       INTEGER garanhao(npopinic),femea(npopinic),randint
       INTEGER vars(npopinic*(kmax-kmin+1)*kmax)
       INTEGER bestvar((kmax-kmin+1)*kmax),kabort
       LOGICAL setk(300),improvement,printfile
       LOGICAL pop(npopinic,300),novager(npopinic,300),rep(300)
       LOGICAL filho(npopinic,300),igual,melhorar
       DOUBLE PRECISION s(300,300),sq(300,300),svector(p*p)
       DOUBLE PRECISION vactual,vcorrente
       DOUBLE PRECISION valores((kmax-kmin+1)*npopinic)
       DOUBLE PRECISION valuepop(npopinic+1),valueger(0:npopinic)
       DOUBLE PRECISION critvalue,amax,aux,somavalue,acumulado,vantes
       DOUBLE PRECISION vmaximo,valuenovager(npopinic)
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

      kabort = kmax+1

      if (printfile) then
        critname(1)='RM'
        critname(2)='RV'
        critname(3)='GCD'
        if (criterio.eq.1) then
          open(3,file='rmgen.txt')
        end if      
        if (criterio.eq.2) then
          open(3,file='rvgen.txt')
        end if      
        if (criterio.eq.3) then
          open(3,file='gcdgen.txt')
        end if
        write(3,*) '        Genetic Algorithm '
        write(3,*) '         Criterion:  ',critname(criterio)
        write(3,*) 'Number of generations requested: ',nger
        write(3,*) ' '
      end if
  
!*********************************
! The loop which is to be repeated for each cardinality of subsets, 
! between kmin e kmax starts here and ends together with the subroutine.
         
      valueger(0)=dfloat(2)
      ordemger(0)=0
      npop=npopinic
      do k=kmin,kmax
         if (printfile) then
           write(3,32) ' CARDINALITY   k=',k
 32        format(' ',/,15x,a20,1X,i3,/)  
         end if
         if (criterio.eq.3) then
             if (.not.esp) then
               nqsi=k
               do m=1,nqsi
                 qsi(m)=m
               end do
             end if
         end if

! Generating the initial population
! Remark: vector rep tests whether all elements of {1,2,...,p} are   
! represented in the initial population, in which case rep(i)=.true. 
! for i=1,...,p
         do i=1,p
            rep(i)=.false.
         end do

         do kpop=1,npopinic
           call randsk1(iseed,p-ndentro,k-ndentro,setk)
           if(ndentro>0) call dcorrigesk(ndentro,dentro,p,setk)
             do i=1,p
               novager(kpop,i)=setk(i)
	       rep(i)=rep(i).or.setk(i)
             end do
         end do

! tests if all elements are represented

         i=1
         rep(p+1)=.false.
         do while(rep(i))
            i=i+1
         end do
         if(i<=p) then 
           write(*,*) '   WARNING:  Not all variables'
           write(*,*) '   are represented in the initial population'
           write(*,*) '   for cardinality k=', k
           write(*,*) '   '
         end if

! Determines the value of the criterion for each member of the initial
! population and ranks them by decreasing order. Vector ordemger will 
! register that ranking. Thus, ordemger(1) will contain the number of the 
! subset in the population with the largest criterion value; ordemger(2)
! will contain the number of the subset with the second largest criterion
! values, and so on.

         somavalue=0
         vmaximo=0
         do kpop=1,npop
           do j=1,p
              setk(j)=novager(kpop,j)
           end do
           if (criterio.eq.1) then
               amax = dobjrm(k,setk,p,s,sq)
               critvalue = dsqrt(amax/tracos)
           end if
           if (criterio.eq.2) then
               amax = dobjrv(k,setk,p,s,sq)
               critvalue = dsqrt(amax/tracosq)
           end if
           if (criterio.eq.3) then
               amax = dobjgcd(nqsi,qsi,valp,vecp,k,setk,p,s,fica)
               critvalue = amax/dsqrt(dfloat(nqsi*k))
           end if

           valueger(kpop) = critvalue

           somavalue=somavalue+valueger(kpop)
           aux=valueger(kpop)
           j=kpop-1
           do while(aux > valueger(ordemger(j)))
             ordemger(j+1) = ordemger(j)
             j=j-1
           end do
           ordemger(j+1) = kpop
         end do

! Orders the initial population by decreasing values of the criterion.

         do kpop=1,npop
           do i=1,p
              pop(kpop,i)=novager(ordemger(kpop),i)
           end do
           valuepop(kpop)=valueger(ordemger(kpop))
         end do


!**************************************************************************
! Here begins a lengthy loop that will repeat the following procedures
! for each generation: arrange marriages, mate and produce offspring, 
! and select, from among the parents 
! and the offspring, the survivors that will make up a new generation.

         do kger=1,nger

!  Select the individuals to make up the couples. The rules are:
!  1. Select npop/2 males ("garanhao"), each male being randomly selected with probability
!      distribution proportional to its criterion value, i.e., 
!      prob(kpop selected)=valuepop(kpop)/sum of valuespop of the population 
!      elements (k-subsets).

           ngaranhao=npop/2
           do kgaranhao=1,ngaranhao
             call slaruv(iseed,1,ranaux)
             r=ranaux*somavalue
             kpop=1
             acumulado=0
             do while(valuepop(kpop)+acumulado<r)
               acumulado=acumulado+valuepop(kpop)
               kpop=kpop+1
             end do
             garanhao(kgaranhao)=kpop
           end do
  
! 2.  For each male, a female ("femea") will be randomly (uniformly) chosen from among the population
!      elements (k-subsets) which differ from the male in at least two genes (variables). This
!      restriction seeks to avoid consanguinity. Note that if the intersection of the parents is of
!      cardinality k-1, offspring (who must inherent the common part of their parents) would
!      necessarily be identical to either the mother or the father.

           nfemea=npop/2
           kfemea=1
           ntalvez=0
           do while(kfemea<=nfemea)
             talvez=randint(iseed,1,npop)
             i=1
             ngenesdistintos=0
             do while (i<=p .and. ngenesdistintos<=1)
	       if (pop(garanhao(kfemea),i).and. .not.pop(talvez,i)) then
                 ngenesdistintos=ngenesdistintos+1
               end if
               i=i+1
             end do
             if (ngenesdistintos>=2) then 
	       femea(kfemea)=talvez
               kfemea=kfemea+1
	       ntalvez=0
	      else
	       ntalvez=ntalvez+1
	       if(ntalvez>10000) then
		  write(*,*) ' '
                  	  write(*,*) ' For cardinality k=', k
                       write(*,*) ' Generation number', kger,' does not'
		  write(*,*) ' have sufficient genetic diversity'
		  write(*,*) ' for acceptable levels of consanguinity' 
		  write(*,*) ' (couples differing by at least 2 genes).' 
		  write(*,*) ' '
		  write(*,*) ' This situation may be overcome by '
                   write(*,*) ' reducing the maximum acceptable number '
                                     write(*,*) ' of clones (nclone).' 
		  write(*,*) ' '
		  write(*,*) ' '
                       write(*,*) ' Best criterion value found so far:',
     +                    valuepop(1)

                  write(*,*) ' ' 
                  kabort = k
                  if (printfile) close(3)
                  return
	       end if   
              end if
            end do

! Reproduction.
! The i-th couple, garanhao(i), femea(i) (with i=1,...,npop/2) will produce a child.
! The child inherits the (possibly empty) common part of both parents. The remaining
! genes (variables in the k-subset) is selected from among (genes unique to the mother) U
! (genes unique to the father). If the cardinality of the mother's unique genes (= cardinality
! of the father's unique genes) is t, k-t genes are randomly (uniformly) selected from among
! the symmetric difference, with the further restriction that not ALL genes will come from a 
! single parent (to avoid clones).
! If the offspring is an identical clone of "nclone" other existing individuals, it will be discarded
! and replaced with a (uniformly) randomly generated k-subset (an adopted child).

            nfilho=ngaranhao
            do kfilho=1,nfilho
              npainaomae=0
              nmaenaopai=0
              do i=1,p
                filho(kfilho,i)=.false.
	       if(pop(garanhao(kfilho),i).and.
     +             .not.pop(femea(kfilho),i)) then
	         npainaomae=npainaomae+1
	         painaomae(npainaomae)=i
	       end if
	       if(.not.pop(garanhao(kfilho),i).and.
     +             pop(femea(kfilho),i)) then
	         nmaenaopai=nmaenaopai+1
	         maenaopai(nmaenaopai)=i
	       end if
	       if(pop(garanhao(kfilho),i).and.
     +             pop(femea(kfilho),i)) then
	           filho(kfilho,i)=.true.
	       end if
              end do
              kpainaomae=randint(iseed,1,npainaomae-1)
              call randsk1(iseed,npainaomae,kpainaomae,setk)
              do i=1,npainaomae
                filho(kfilho,painaomae(i))=setk(i)
              end do
              kmaenaopai=npainaomae-kpainaomae
              call randsk1(iseed,nmaenaopai,kmaenaopai,setk)
              do i=1,nmaenaopai
                filho(kfilho,maenaopai(i))=setk(i)
              end do

! Determines the (nvezes) number of times which this genetic makeup exists in the population
! (including newborn children). If nvezes > nclone, the child is rejected and a new child of unknown
! parents is adopted.

              nvezes=0
	      kpop=1
	      do while(nvezes<=nclone .and. kpop<=npop .and. 
     +            npop-kpop+kfilho>nclone-nvezes)
	        j=1
	        igual=.true.
	        do while(j<=p .and. igual)
	          if(filho(kfilho,j).NEQV.pop(kpop,j)) then
		     igual=.false.
                   else
                     j=j+1
	          end if
	        end do
	        if(igual) nvezes=nvezes+1
	        kpop=kpop+1
	      end do
	      kfilhoant=1
	      do while(nvezes<=nclone .and. kfilhoant<kfilho .and. 
     +            kfilho-kfilhoant>=nclone-nvezes)
	        j=1
	        igual=.true.
	        do while(j<=p .and. igual)
	          if (filho(kfilho,j).NEQV.filho(kfilhoant,j)) then
		      igual=.false.
                   else
                      j=j+1
	          end if
	        end do
	        if(igual) nvezes=nvezes+1
	        kfilhoant=kfilhoant+1
	      end do
              if(nvezes>nclone) then
              call randsk1(iseed,p-ndentro,k-ndentro,setk)
              if(ndentro>0) call dcorrigesk(ndentro,dentro,p,setk)
              do j=1,p
                filho(kfilho,j)=setk(j)
              end do
             end if
            end do

! Calculates each child's value of the criterion, and ranks the children accordingly (descending
! order). Thus, ordemger(1) will be the child with the largest value of the criterion, ordemger(2)
! the child with the second largest criterion value, etc.

            do kfilho=1,nfilho
              do j=1,p
	        setk(j)=filho(kfilho,j)
              end do
              if (criterio.eq.1) then
                 amax = dobjrm(k,setk,p,s,sq)
                 critvalue = dsqrt(amax/tracos)
              end if
              if (criterio.eq.2) then
                 amax = dobjrv(k,setk,p,s,sq)
                 critvalue = dsqrt(amax/tracosq)
              end if
              if (criterio.eq.3) then
                 amax = dobjgcd(nqsi,qsi,valp,vecp,k,setk,p,s,fica)
                 critvalue = amax/dsqrt(dfloat(nqsi*k))
              end if
              valueger(kfilho)=critvalue

! If the logical variable "melhorar" is .true., the newborn child will be genetically improved (i.e.,
! subject to the restricted local improvement algorithm). 
! WARNING: This option may tend to generate large numbers of clones. Default is .false.

              if(melhorar) then
                vantes=amax
                call dmelhoramentogen(criterio,p,setk,amax,ndentro,
     +              dentro,k,s,sq,nqsi,qsi,valp,vecp,fica)
                if(amax>vantes) then
                  if (criterio.eq.1) then
                    critvalue = dsqrt(amax/tracos)
                  end if
                  if (criterio.eq.2) then
                    critvalue = dsqrt(amax/tracosq)
                  end if
                  if (criterio.eq.3) then
                    critvalue = amax/dsqrt(dfloat(nqsi*k))
                  end if
                  valueger(kfilho)=critvalue
                  do j=1,p
                    filho(kfilho,j)=setk(j)
                  end do
                end if
              end if
 
              aux=valueger(kfilho)
              j=kfilho-1
              do while(aux>valueger(ordemger(j)))
                 ordemger(j+1)=ordemger(j)
	         j=j-1
              end do
              ordemger(j+1)=kfilho
            end do

! Selection of a new generation.
! The new generation is made up of the npop best (highest criterion value) elements
! from among the parents and their offspring. novager(1,i), (i=1,p) is the best element 
! inthe new generation, novager(2,i) (i=1,p) is the second best, and so on.

            kfilho=1
            kpop=1
            knovager=1
            do while(knovager<=npop)
              if(kfilho>nfilho) then
                do while(knovager<=npop)
	          do i=1,p
	            novager(knovager,i)=pop(kpop,i)
	          end do
	          valuenovager(knovager)=valuepop(kpop)
                  kpop=kpop+1
	          knovager=knovager+1
	        end do
               else
	        if(valueger(ordemger(kfilho))>valuepop(kpop)) then
	          do i=1,p
	            novager(knovager,i)=filho(ordemger(kfilho),i)
	          end do
	          valuenovager(knovager)=valueger(ordemger(kfilho))
	          kfilho=kfilho+1
	          knovager=knovager+1
	         else
	          do i=1,p
	           novager(knovager,i)=pop(kpop,i)
	          end do
	          valuenovager(knovager)=valuepop(kpop)
                  kpop=kpop+1
	          knovager=knovager+1
	        end if
              end if
            end do
   
! The survivors (i.e., the new population) make up the next generation.
            
            do kpop=1,npop
              do i=1,p
                pop(kpop,i)=novager(kpop,i)
              end do
              valuepop(kpop)=valuenovager(kpop)
            end do

            if(valuepop(1)>vmaximo) then
              vmaximo=valuepop(1)
              if (printfile) then
                write(3,20) ' In generation',kger,'a new maximum =',
     +               valuepop(1)
 20             format(' ',a15,i4,1x,a9,f12.7)
              end if
            end if
         end do 

         if (printfile) then
           write(3,24) ' Best solution found after',nger,' generations:'
 24        format(' ',/,a28,1x,i4,a13)
         end if

! auxiliary vector for a one-line print-out
         naux=0
         do j=1,p
           if(pop(1,j)) then
             naux=naux+1
             auxw(naux)= fica(j)
           end if
         end do

         critvalue=valuepop(1)

         if (printfile) then 
           write(3,*) ' ',(auxw(j), j=1,naux)
           write(3,33) ' Best value of the ',critname(criterio),
     +      ' criterion  =', critvalue
 33      format(' ',a20,1x,a3,a13,f10.7)
         end if

!  preparing the output which will be used in R.
          do kpop=1,npop
           jjaux=0
           do i=1,p
             if (pop(kpop,i)) then
               jjaux=jjaux+1
               vars(npop*kmax*(k-kmin)+npop*(jjaux-1)+kpop) = fica(i)
             end if
           end do  
           valores((k-kmin)*npop+kpop) = valuepop(kpop) 
          end do

! Optional further improvement of the best solution
         if (improvement) then
           if (criterio.eq.1) vactual=valuepop(1)**2*tracos
           if (criterio.eq.2) vactual=valuepop(1)**2*tracosq
           if (criterio.eq.3) vactual=valuepop(1)*dsqrt(dfloat(nqsi*k))
           do j=1,p
            setk(j)=pop(1,j)
           end do
           vcorrente = vactual

           call dmelhoramentogen(criterio,p,setk,vactual,ndentro,dentro,
     +      k,s,sq,nqsi,qsi,valp,vecp,fica)


           if (vactual > vcorrente) then 
                   naux=0
                   do j=1,p
                     if(setk(j)) then
                       naux=naux+1
                       auxw(naux)= fica(j)
                     end if
                   end do
                   if (criterio.eq.1) then
                       critvalue=dsqrt(vactual/tracos)
                   end if
                   if (criterio.eq.2) then
                       critvalue=dsqrt(vactual/tracosq)
                   end if
                   if (criterio.eq.3) then
                       critvalue = vactual/dsqrt(dfloat(nqsi*k))
                   end if
                  if (printfile) then
                     write(3,25) ' Improved ',critname(criterio), 
     +              'value =',critvalue, '  Improved subset:',
     +              (auxw(j),j=1,naux)
 25                format(' ',/,a12,a3,a8,1X,f10.7,/,1X,a18,/,1x,18(i4)) 
                  end if
            else      
                  if (printfile) then
                     write(3,*) 'The solution was not improved'
                  end if
            end if      
         end if

! End of eventual improvement of the best solution

         bestval(k-kmin+1)=critvalue
         do j=1,k
            bestvar(kmax*(k-kmin)+j)=auxw(j)
         end do
         end do 

! End of k-loop

         if (printfile) close(3)
         return
       end



