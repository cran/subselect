* SUBPROGRAMS:
* Included in the R/S package "subselect", available on CRAN.
********************************************************************
       subroutine newinicializar(criterio,p,s,svector,sq,nfora,
     +  fora,ndentro,dentro,fica,tracos,tracosq,vecp,poriginal,vecvecp,
     +   h,hvector,rh)
**********************************************************************
* performs some initial chores: converts svector into a matrix s; computes
* the square of that matrix, sq, and the traces of both s (for RM) and sq 
* (for RV); converts the eigenvectors of s (for GCD), calculated within R 
* and passed on as vector vecvecp, into a matrix. All this is done taking 
* into account the variables that are to be 
* forcibly excluded (vector fora) or included (vector dentro).
*
*
* INPUT:
*
* criterio - integer variable indicating the criterion for the quality of 
*            the k-subset.
*    p     - integer variable, giving the number of original variables.
* svector  - double precision vector, giving the full covariance (or 
*            correlation) matrix of the p variables (given as a vector
*            for convenience in passing from R to Fortran). 
*  nfora   - integer variable indicating the number of original variables 
*            that are to be forcefully excluded from the k-subset.
*   fora   - integer vector, indicating the variable numbers associated
*            with the nfora variables that are to be excluded by force.
* ndentro  - integer variable indicating the number of original variables 
*            that are to be forcefully included in the k-subset.
*  dentro  - integer vector, indicating the variable numbers associated
*            with the ndentro variables that are to be included by force.
* vecvecp  - double precision vector, giving the eigenvectors of the 
*            covariance (or correlation) matrix of the p variables (given 
*            as a vector for convenience in passing from R to Fortran).
* hvector  - double precision vector, giving the full effect descrption 
*            matrix (H) of the p variables (given as a vector for
*            convenience in passing from R to Fortran). 
*   rh     - H matrix rank (expected)
*
* OUTPUT:
*
*   p      - integer variable indicating the number of admissible variables,
*            i.e., the original number of variables p minus nfora.
*   s      - double precision array (pxp) giving the covariance (or 
*            correlation) matrix of the admissible variables.
*   sq     - double precision array (pxp) giving the square of matrix s.
*  fica    - integer vector giving the original variable numbers of the 
*            admissible variables.
* tracos   - double precision variable giving the trace of matrix s.
* tracosq  - double precision variable giving the trace of matrix sq.
*  vecp    - double precision array (pxp) giving the eigenvectors of s (in
*            the same column-order as the eigenvalues).
* poriginal - integer variable, indicating the original number of variables
*             (poriginal = p + nfora)
*   h       - double precision array (pxp) giving the full effect descrption 
*            matrix (H)  
*

* general declarations
      integer p,poriginal,nfora,ndentro,iaux,nfica,criterio,rh
      integer fora(0:nfora),fica(0:p),dentro(0:ndentro)
*     integer dm
      double precision s(p,p),sq(p,p),svector(p*p),h(p,p)
	double precision vecvecp(p*p),hvector(p*p)
* declaration only for the RM criterion
      double precision tracos
* declarations only for the RV criterion
      double precision tracosq
* declarations only for the GCD criterion
      double precision vecp(p,p)
*      double precision valp(p)
* integer jaux

      external dprodmat



* conversion of svector into a matrix (s); if GCD, initialization of the vecp
* matrix which will hold the eigenvectors of s (computed in R and passed 
* in as vector vecvecp)
       do j=1,p
         DO i= 1,p
           s(i,j) = svector(i+(j-1)*p)
           if (criterio.eq.3) then 
              vecp(i,j) = vecvecp(i+(j-1)*p)
           end if 
         END DO
      end do

       If (rh .GT. 0) then
* conversion of hvector into a matrix (h) 
       do j=1,p
         DO i= 1,p
           h(i,j) = hvector(i+(j-1)*p)
          END DO
        end do
        end if


* Defines matrix sq := s x s. Computes the traces of s (for RM) and 
* sq (for RV), and the eigendecomposition of s (for GCD)


      call dprodmat(p,p,s,p,s,sq)


      if (criterio .eq. 1) then 
         tracos=0
         DO j=1,p
           tracos=tracos+s(j,j)
         END DO
      end if
      if (criterio .eq. 2) then 
         tracosq=0
         DO j=1,p
           tracosq=tracosq+sq(j,j)
         END DO
      end if


* Re-defines the matrices s, sq and h (if r>0), excluding the rows/columns that the 
* user has requested be excluded from the final solution. 
* p is redefined as p-nfora, where nfora is the number of variables that 
* are to be excluded. The original variable-number for each of the surviving
* variables is preserved in vector fica(j), j=1,p-nfora.
* Warning: the original matrices s and sq are destroyed.

      poriginal=p
      p=p-nfora
      fora(0)=0
      nfica=0
      fica(0)=0
      do j=0,nfora-1
         do i=fora(j)+1,fora(j+1)-1
            nfica=nfica+1
            fica(nfica)=i
         end do
      end do
      iaux=fora(nfora)
      do j=iaux+1,poriginal
         nfica=nfica+1
         fica(nfica)=j
      end do

      if(nfora .GT. 0) then
         do i=1,p-1
            do j=i+1,p
              s(i,j)=s(fica(i),fica(j))
              s(j,i)=s(i,j)
              sq(i,j)=sq(fica(i),fica(j))
              sq(j,i)=sq(i,j)
              if(rh .GT. 0) then
		    h(i,j)=h(fica(i),fica(j))
              h(j,i)=h(i,j)
		    end if
		  end do
            s(i,i)=s(fica(i),fica(i))
            sq(i,i)=sq(fica(i),fica(i))
		    if(rh .GT. 0) then
		    h(i,i)=h(fica(i),fica(i))
              end if
         end do
         s(p,p)=s(fica(p),fica(p))
         sq(p,p)=sq(fica(p),fica(p))
          if(rh .GT. 0) then
		    h(p,p)=h(fica(p),fica(p))
          end if 

      end if

* The vector of variables that are forced to be included in the solution 
* is redefined, taking into account that nfora variables (may) have been 
* excluded.

      if((ndentro.GT.0) .and. (nfora.GT.0)) then
         i=1
         j=1
         do while(j.LE.ndentro)
           if(dentro(j).EQ.fica(i)) then
             dentro(j)=i
             j=j+1
           end if
           i=i+1
         end do
      end if
      return
      end
**********************************************************************

*********************************************************************
      subroutine randsk1(n,k,sk)
**********************************************************************
* generates a subset of k random integers from the set {1,2,...,n}. 
* Warning: uses function randint.

*  INPUT: 
*    iseed
*    n   - integer variable, largest integer in the set.
*    k   - integer variable, giving the number of random integers that
*                are to be generated
* 
*  OUTPUT: 
*   sk   - logical vector defining the generated subset 
*               (sk(i)=.true. iff i belongs to sk, i=1,...,n)

       integer randint,pp(n)
       logical sk(n)



       do i=1,n
        sk(i)=.false.
        pp(i)=i
       end do
       do i=1,k
        nalt=randint(i,n)
        sk(pp(nalt))=.true.
        pp(nalt)=pp(i)
       end do
       return 
      end
***********************************************************************
************************************************************************
      integer function randint(esq,dir)
**********************************************************************
* randomly generates a number in the interval [esq,dir].
*
* WARNING: Uses a call to R's random number generator (as described
*  in pages 48 and 50-51 of the "Writing R Extensions" manual, version
*  1.5).
*
*INPUT: 
       
*        esq  - integer variable, giving lower bound of interval.
*        dir  - integer variable, giving upper bound of interval.

       external unifrnd, rndstart, rndend
       
       double precision ran
*       real ran
       integer esq,dir


       call rndstart()
       ran = unifrnd()
       randint=esq+int(ran*(dir-esq+1))
       call rndend()
       return 
      end
*********************************************************************



*********************************************************************
      subroutine dprodmat(n,m,a,r,b,prod)
**********************************************************************
* calculates the matrix product of two matrices
*
* INPUT : n - integer, number of rows of matrix a
*         m - integer, number of columns of matrix a (= number of rows
*                matrix b)
* 
*         a - double precision array (nxm)
*         r - integer, number of rows of matrix b
*         b - double precision array (mxr)
* OUTPUT: 
*   prod - double precision array (nxr), the product a x b.

       INTEGER r
       DOUBLE PRECISION a(n,m),b(m,r),prod(n,r),soma

       DO i=1,n
         DO j=1,r
           soma=0 
        	DO l=1,m
	          soma=soma+a(i,l)*b(l,j)
	        END DO
	        prod(i,j)=soma
         END DO
       END DO
       RETURN
      END
**********************************************************************

**********************************************************************
       subroutine dcorrigesk(ndentro,dentro,p,setk,poriginal)
**********************************************************************
* subroutine randsk1, randomly selects (with uniform 
* distribution) a subset (setk) with cardinality k-ndentro  of the set 
* {1,2,...,p-ndentro}. This routine takes that subset and merges it with 
* the ndentro variables indicated in vector dentro to produce a subset of
* cardinality k=|setk|+ndentro, whose variable numbers are in accordance with 
* the original variable numbers. The vector setk is re-defined accordingly.
* Thus, the two subroutines randsk1and dcorrigesk randomly select a subset
* of cardinality k from the set {1,2,...,p}, ensuring that the ndentro 
* variables specified by dentro are included.
*
* INPUT: 
*
* ndentro - integer variable giving the size of vector dentro.
* dentro  - integer vector giving the variable numbers of variables that are
*           to be forcefully included in the k-subset.
*   p     - number of elements in the full set (excluding 'exclude' variables).
*  setk   - logical vector with setk(j)=.true. iff j belongs to the  
*           subset of {1,2,...,p-ndentro}
* poriginal - number of variables originally (equal to 'p', except when
*              the 'exclude' option has been used
* 
* OUTPUT: 
* 
*  setk   - logical vector with setk(j)=.true. iff j belongs to the  
*           subset of {1,2,...,p}

       integer dentro(0:ndentro),p,poriginal
       logical setk(poriginal)

       dentro(0)=0
       naux=ndentro
       i=p-ndentro
       j=p
       do while(j.GE.1)
         if(dentro(naux).EQ.j) then
           setk(j)=.true.
	   naux=naux-1
	 else
	   setk(j)=setk(i)
	   i=i-1
         end if
         j=j-1
       end do
       return
      end
**********************************************************************

************************************************************************
      subroutine dannealing(criterio,p,k,s,sq,setk,vactual,ndentro,
     +  dentro,niter,beta,temp,coolfreq,tracos,tracosq,nqsi,qsi,
     +  valp,vecp,fica,poriginal,h,rh)
**********************************************************************
* 
* Applies simulated annealing to a subset setk of the original variables. 
*
*
* INPUT : 
*
*criterio - integer variable indicating the criterion of subset quality
*            requested.
*   p     - integer variable giving the number of original variables (minus
*               ndentro).
*   k     - integer variable, cardinality of subset.
*   s     - double precision matrix of covariances of all p variables.
*   sq    - double precision matrix, the square of s. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*ndentro  - integer variable, giving the number of variables that are to
*          be forcefully included in the output subset.
* dentro  - integer vector, indicating the numbers associated with the 
*          ndentro variables that are to be forced into the subsets.
* niter   - integer variable, number of simulated annealing iterations 
*           requested for each solution.
*  beta   - double precision variable,(between 0 and 1) indicating the
*           geometric cooling factor for simulated annealing.
* temp   - double precision variable, indicating the initial temperature
*           to e used in the simulated annealing algorithm.
* coolfreq - integer indicating the cooling frequency, i.e., the number 
*           of iterations that are to be carried out for each fixed
*           temperature 
* tracos  - double precision variable, the trace of the covariance (or
*           correlation) matrix s. (For RM criterion only).
* tracosq - double precision variable, the trace of the matrix sq. (For RV
*           criterion only).
*  nqsi   - integer variable, indicating the dimension of the Principal 
*           Subspace with which the subspace spanned by the k-subset is
*            to be compared. (For GCD criterion only).
*  qsi    - integer vector, giving the ranks of the Principal Components 
*           which span the Principal Subspace described above (see nqsi).
*           (For GCD criterion only).
*  valp   - double precision vector giving the eigenvalues of s, in 
*            decreasing order. (For GCD criterion only).
*  vecp   - double precision array (pxp) giving the eigenvectors of s, in
*            the same column-order as the values of valp. (For GCD only).
*  fica   - integer vector giving the original variable numbers of the 
*            admissible variables.
* poriginal - 
*
*   h     - double precision matrix giving the full effect descrption 
*            matrix (H)  
*
*   rh     - integer variable giving the expected rank of H matrix.


*
* OUTPUT: 
*
* setk    - logical vector indicating the best subset produced by simulated 
*           annealing after the niter iterations. setk(j)=.true. iff variable 
*           j belongs to the subset. 
*vactual  - double precision variable giving the best value of the criterion 
*           produced throughout the niter iterations.


        external unifrnd, rndstart, rndend
        
* general declarations

* character*3 aceita
        double precision r
        integer coolfreq, criterio, p, poriginal
        integer dentro(0:ndentro),que(p),cons(p),randint
        logical setk(poriginal),setkmax(poriginal),auxlog(p)
        double precision s(poriginal,poriginal),h(poriginal,poriginal)
        double precision sq(poriginal,poriginal),beta,citer,temp
        double precision vactual,vmax,vtroca,dir
* double precision crittroca,critactual
* declarations specific to the RM criterion
        double precision tracos,dobjrm
* declarations specific to the RV criterion
        double precision tracosq,dobjrv
* declarartions specific to the GCD criterion
        integer qsi(p),fica(0:p),nqsi
        double precision valp(poriginal)
        double precision vecp(poriginal,poriginal),dobjgcd
* declarations for the tau_2, xi_2 criteria
        double precision dobjtau2,dobjxi2,dobjzeta2,dobjccr12
        integer rh

   

* initializing 
        citer = temp
        vtroca = 0.0
        dir = 0.0

        do j=1,p
         auxlog(j)=.true.
        end do
        do j=1,ndentro
         auxlog(dentro(j))=.false.
        end do

        ncons=0
        nque=0
        do j=1,p
           if(.not.setk(j)) then
              nque=nque+1
	      que(nque)=j
	   else
	      if(auxlog(j)) then
	        ncons=ncons+1
	        cons(ncons)=j
	      end if
           end if
        end do

        vmax=vactual
        do j=1,p
          setkmax(j)=setk(j)
        end do

* Leftovers of previous tests, if anyone should ever be interested in
* seeing a HUGE file of subset values..., all subsequent commented code 
* beginning with "if(printall)" used 
* to print to file *ALL* solutions considered by simulated annealing 
*        if (printall) then
*          open(16,file='safull.out')
*          write(16,*) '    iter','  citer','    vtroca','     vactual',  
*     +     'naleat [0,1]','    moeda ao ar','   aceita'
*        end if

* more initializing
        iter=0      
        nigual=0
        inaoconv=0
        do while(iter.LT.niter)
         
 
* generates a (uniform) random integer, kque, in {1,2,...,nque}.
* jentra, the variable which is being considered for inclusion in the subset
* setk, is the element in position kque of the vector que 

         kque=randint(1,nque)
         jentra=que(kque)
  
* generates a (uniform) random integer, kcons, in {1,2,...,kcons}.
* jsai, the variable which is being considered for exclusion from the subset
* setk, is the element in position kcons of the vector cons, i.e., an 
* element of setk not belonging to vector dentro

         kcons=randint(1,ncons)
         jsai=cons(kcons)

* calculates the value of the criterion for solution setk\{jsai}U{jentra}

         setk(jsai)=.false.
         setk(jentra)=.true.
         if (criterio.eq.1) then
           vtroca=dobjrm(k,setk,p,s,sq,poriginal)
         end if
         if (criterio.eq.2) then
           vtroca=dobjrv(k,setk,p,s,sq,poriginal)
         end if
         if (criterio.eq.3) then
           vtroca=dobjgcd(nqsi,qsi,valp,vecp,k,setk,p,s,fica,poriginal)
         end if
         if (criterio.eq.4) then
           vtroca=dobjtau2(k,setk,p,s,poriginal,h,rh)
         end if
	 if (criterio.eq.5) then
           vtroca=dobjxi2(k,setk,p,s,poriginal,h,rh)
         end if
	 if (criterio.eq.6) then
           vtroca=dobjzeta2(k,setk,p,s,poriginal,h,rh)
         end if
	 if (criterio.eq.7) then
           vtroca=dobjccr12(k,setk,p,s,poriginal,h,rh)
         end if
* if vtroca>vmax, updates vmax and setkmax, which hold the best value and 
* subset found so far.

         if(vtroca.GT.vmax) then
          vmax=vtroca
          do j=1,p
	    setkmax(j)=setk(j)
	  end do
         end if

* The replacement is carried out if vtroca>=vactual or, if not, with 
* probability given by the simulated annealing algorithm, and computed
* by R's Random Number Generator


         call rndstart()
         r = unifrnd()
         call rndend()

         if(vtroca-vactual.GE.0) then 
          dir=10000.0D0
         else
          if (criterio.eq.1) then
           dir=exp((dsqrt(vtroca)-dsqrt(vactual))/(dsqrt(tracos)*citer))
          end if
          if (criterio.eq.2) then
          dir=exp((dsqrt(vtroca)-dsqrt(vactual))/(dsqrt(tracosq)*citer))
          end if
          if (criterio.eq.3) then
          dir=exp(((vtroca/dsqrt(DBLE(nqsi*k)))-
     +         (vactual/dsqrt(DBLE(nqsi*k))))/citer)
          end if
	    if (criterio.eq.4) then
          dir=exp((vtroca-vactual)/citer)
          end if
		if (criterio.eq.5) then
          dir=exp((vtroca-vactual)/citer)
          end if
		if (criterio.eq.6) then
          dir=exp((vtroca-vactual)/citer)
          end if
		if (criterio.eq.7) then
          dir=exp((vtroca-vactual)/citer)
          end if
         end if
         if(vtroca.GE.vactual.or.(vtroca.LT.vactual.and.r.LT.dir)) then 
* swap

          vactual=vtroca
          que(kque)=jsai
          cons(kcons)=jentra

*           if (printall) aceita='sim'
*            aceita = 'sim'

           nigual=0
        else 
* reject the swap
           setk(jsai)=.true.
           setk(jentra)=.false.

*           if (printall) aceita='nao'
*           aceita = 'nao'

          nigual=nigual+1 
* counts number of consecutive rejections of a swap
        end if

* Hacking the temperature:
* updates the temperature citer every coolfreq iterations

        if(mod(iter,coolfreq).EQ.0) citer=citer*(1-beta)


* when there are more than p consecutive rejections of a swap and the
* first 20p iterations have not yet been completed, the
* simulated annealing temperature citer will be increased: by a factor 
* of 2 if less than 5p iterations have gone by, a factor of 1.5 if the
* number of iterations is between 5p and 10p, by a factor of 1.1 if
* the current number of iterations is between 10p and 20p. 

       if(nigual.GE.p) then
         if(iter.LE.5*p) then
              citer=citer*2.0D0
	 else
	      if(iter.LE.10*p) then
	        citer=citer*1.5D0
	      else
	        if(iter.LE.20*p) then
	            citer=citer*1.1D0
	        else
		    inaoconv=inaoconv+1
	        end if
	       end if
          end if
          nigual=0
        end if

* inaoconv counts the number of times in which the currently best solution 
* has remained unchanged after the first 20p iterations.

* inaoconv counts the number of consecutive iterations without changes in
* the current solution AFTER THE FIRST 20p ITERATIONS. If this becomes too
* large (the value chosen here was 500), the algorithm stops (a further change
* is unlikely, given the low temperature after so many iterations).

         if (inaoconv.GT.500) iter=niter
         iter=iter+1
          
*         if (printall) then 
*           if (criterio.eq.1) then                 
*             crittroca=dsqrt(vtroca/tracos)        
*             critactual=dsqrt(vactual/tracos)      
*           endif                                   
*           if (criterio.eq.2) then                 
*             crittroca=dsqrt(vtroca/tracosq)      
*             critactual=dsqrt(vactual/tracosq)     
*           endif                                    
*           if (criterio.eq.3) then                 
*             crittroca=vtroca/dsqrt(DBLE(nqsi*k)) 
*             critactual=vactual/dsqrt(DBLE(nqsi*k))
*           endif                                     
*            write(16,101) iter,citer,crittroca,      
*     +       critactual,r,dir,aceita
*101         format(' ',1X,i6,1X,f7.4,2x,f10.7,1X,f10.7,7x,f6.4,
*     +       9x,f16.14,7x,a3)                        
*         end if


      end do
 
* defines the algorithm's solution, vactual, as the best solution obtained
* at any stage of the niter iterations. 

       vactual=vmax
       do i=1,p
          setk(i)=setkmax(i)
       end do


       return
      end



**********************************************************************

**********************************************************************
       subroutine dmelhoramentogen(criterio,p,setk,vactual,ndentro,
     +  dentro,k,s,sq,nqsi,qsi,valp,vecp,fica,poriginal,h,rh)
**********************************************************************
* This subroutine (which is only called by the "improve" subroutine, 
* but also by the other two main routines - "genetic" and "anneal" -  if 
* the logical variable "improvement" is set to .true.) seeks to improve an 
* initial k-subset by a modified local search algorithm, the details of 
* which are as follows. The variables not belonging to this initial subset 
* are placed in a queue ("que"). This subroutine explores the possibility 
* of replacing a variable in the subset with a variable from the queue. 
* More precisely, a variable j is selected and removed from the queue. Each 
* variable i in the subset is, in turn, replaced by variable j and the 
* resulting values of the criterion are computed. If the best of these k 
* new criterion values exceeds the subset's original criterion value, the 
* current solution is updated accordingly. In this case, the variable which
* leaves the subset is added to the queue, but only if it has not previously 
* been in the queue (i.e., no variable can enter the queue twice). The 
* algorithm proceeds until the queue is emptied. 
*
* INPUT: 
*criterio - integer variable indicating the criterion of subset quality
*            requested.
*   p     - integer variable, number of original variables. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on exit. 
*vactual  - double precision variable giving the best value of the criterion 
*           for the initial solution. Changed on exit.
*ndentro  - integer variable, giving the number of variables that are to
*          be forcefully included in the output subset.
* dentro  - integer vector, indicating the numbers associated with the 
*          ndentro variables that are to be forced into the subsets.
*   k     - integer, cardinality of subset setk.
*   s     - double precision matrix of covariances of all p variables.
*   sq    - double precision matrix, the square of s.
*  nqsi   - integer variable, indicating the dimension of the Principal 
*           Subspace with which the subspace spanned by the k-subset is
*            to be compared. (For GCD criterion only).
*  qsi    - integer vector, giving the ranks of the Principal Components 
*           which span the Principal Subspace described above (see nqsi).
*           (For GCD criterion only).
*  valp   - double precision vector giving the eigenvalues of s, in 
*            decreasing order. (For GCD criterion only).
*  vecp   - double precision array (pxp) giving the eigenvectors of s, in
*            the same column-order as the values of valp. (For GCD only).
*  fica   - integer vector giving the original variable numbers of the 
*            admissible variables.
* poriginal - 
*
*   h     - double precision matrix giving the full effect descrption 
*            matrix (H)  
*
*   rh     - integer variable giving the expected rank of H matrix.
*
*
* OUTPUT: 
*
* setk    - logical vector indicating the best subset after the modified local
*           search, with setk(j)=.true. iff variable j belongs to the subset. 
*vactual  - double precision variable giving the best value of the criterion 
*           after the modified local search.

* general declarations
      integer p,poriginal,jmax,jconsmax
      integer dentro(0:ndentro),que(p),cons(p),ndentro,k,criterio
      logical setk(poriginal),auxlog(p),esteveque(p)
      double precision s(poriginal,poriginal),h(poriginal,poriginal)
      double precision vactual,vtroca,vtrocamax
* declarations for the RM and RV criteria
      double precision sq(poriginal,poriginal),dobjrm,dobjrv
* declarations for the gcd criteria
        integer qsi(p),fica(0:p),nqsi
        double precision valp(poriginal)
        double precision vecp(poriginal,poriginal),dobjgcd
* declarations for the tau_2, xi_2 criteria
        double precision dobjtau2,dobjxi2,dobjzeta2,dobjccr12
        integer rh


* initializations
      vtroca = 0.0
      jmax = 0 
      jconsmax = 0

      do j=1,p
        auxlog(j) = .true.
      end do
      do j=1,ndentro
        auxlog(dentro(j))=.false.
      end do

      ncons=0
      nque=0
      DO j=1,p
         IF (.not.setk(j)) THEN
           nque=nque+1
	   que(nque)=j
	   esteveque(j)=.true.
         ELSE
	   esteveque(j)=.false.
	   IF (auxlog(j)) then
             ncons=ncons+1
             cons(ncons)=j
	   END IF
         END IF
      END DO

* emptying the queue "que"
       DO WHILE(nque.GT.0 .and. ncons.GT.0)
       jentra=que(nque)
       nque=nque-1
       vtrocamax=0
       DO j=1,ncons
        jsai=cons(j)
	setk(jsai)=.false.
	setk(jentra)=.true.
        if (criterio.eq.1) then
        	vtroca=dobjrm(k,setk,p,s,sq,poriginal)
        end if
        if (criterio.eq.2) then
        	vtroca=dobjrv(k,setk,p,s,sq,poriginal)
        end if
        if (criterio.eq.3) then
                vtroca=dobjgcd(nqsi,qsi,valp,vecp,k,setk,p,s,
     +              fica,poriginal)
        end if
	  if (criterio.eq.4) then
                vtroca=dobjtau2(k,setk,p,s,poriginal,h,rh)
        end if
	  if (criterio.eq.5) then
                vtroca=dobjxi2(k,setk,p,s,poriginal,h,rh)
        end if
	    if (criterio.eq.6) then
                vtroca=dobjzeta2(k,setk,p,s,poriginal,h,rh)
        end if
	  if (criterio.eq.7) then
                vtroca=dobjccr12(k,setk,p,s,poriginal,h,rh)
        end if
	
	IF (vtroca.GT.vtrocamax) THEN
           vtrocamax=vtroca
           jmax=jsai
           jconsmax=j
	END IF
	setk(jsai)=.true.
	setk(jentra)=.false.
       END DO
       IF (vtrocamax.GT.vactual) THEN
        vactual=vtrocamax
        setk(jmax)=.false.
	setk(jentra)=.true.
	cons(jconsmax)=jentra
        IF (.not.esteveque(jmax)) THEN
           nque=nque+1
           que(nque)=que(1)
           que(1)=jmax
           esteveque(jmax)=.true.
	END IF
       END IF
       END DO
       RETURN
       END
**********************************************************************


**********************************************************************
      function dobjrm(k,setk,p,s,sq,poriginal)
**********************************************************************
* This function computes the "variable part" of the RM criterion. By "variable
* part" is meant that which changes from one given k-subset to the other 
* (the final value of RM is only computed in the main subroutine, for 
* the selected solutions, in order to spare some computations in this
* subroutine, which is called extensively). 
*
* WARNING : This function calls the LAPACK routine DPOSV.
*
*  OUTPUT:
*
* The trace of the matrix inv(S_k) x (S**2)_k, where S_k 
* is the submatrix of the covariance (or correlation) matrix for the
* variables in subset setk, (S**2)_k is its equivalent for the square of
* that covariance/correlation matrix, and "inv" stands for matrix inverse.
*
*  INPUT: 
*
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables.
*   sq    - double precision matrix, the square of s.
* poriginal - 


       external dposv

       integer p,poriginal,setint(p)
*      integer dm
       integer i,j,k
       logical setk(poriginal)
       character*1 laux
       double precision s(poriginal,poriginal),sq(poriginal,poriginal)
       double precision sk(k,k),skinput(k,k),dobjrm

* initializations and call to LAPACK routine DPOSV
* matrix skinput will store the submatrix S_k. Matrix s will initially
* be identity, fed to the LAPACK routine DPOSV, from which it will emerge
* as the inverse of S_k. 


       do i=1,p
          setint(i)=i
       end do


       i=0
       do j=1,p
         if(setk(j)) then
           i=i+1
           setint(i)=j
         end if
       end do

       do i=1,k-1
           do j=i+1,k
              skinput(i,j)=s(setint(i),setint(j))
              skinput(j,i)=skinput(i,j)
              sk(i,j)=0
              sk(j,i)=0
           end do
           skinput(i,i)=s(setint(i),setint(i))
           sk(i,i)=1
         end do

          skinput(k,k)=s(setint(k),setint(k))
          sk(k,k)=1

*          dm=300
          laux='L'
          info=0
          call dposv(laux,k,k,skinput,k,sk,k,info)

* computes the trace of inv(S_k) x (S**2)_k

          dobjrm=0
	  do i=1,k
	    do l=1,k
		  dobjrm=dobjrm+sk(i,l)*sq(setint(l),setint(i))
            end do
	  end do
	  return
	end   
**********************************************************************

**********************************************************************
      function dobjrv(k,setk,p,s,sq,poriginal)
**********************************************************************
* This function computes the "variable part" of the RV criterion. By "variable
* part" is meant that which changes from one given k-subset to the other 
* (the final value of RV is only computed in the main subroutine, for 
* the selected solutions, in order to spare some computations in this
* subroutine, which is called extensively). 
*
* WARNING : This function calls the LAPACK routine DPOSV.
* 
* OUTPUT:
* The trace of the matrix (inv(S_k) x (S**2)_k)**2, where S_k 
* is the submatrix of the covariance (or correlation) matrix for the
* variables in subset setk, (S**2)_k is its equivalent for the square of
* that covariance/correlation matrix, and "inv" stands for matrix inverse.
*
*  INPUT: 
*
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables.
*   sq    - double precision matrix, the square of s.
* poriginal -

      external dposv

      integer p,poriginal,setint(p),info
* integer dm
      logical setk(poriginal)
      character*1 laux
      double precision s(poriginal,poriginal),sq(poriginal,poriginal)
      double precision sk(k,k),skinput(k,k),soma1,soma2,soma,dobjrv

* initializations and call to LAPACK routine DPOSV
* matrix skinput will store the submatrix S_k. Matrix s will initially
* be identity, fed to the LAPACK routine DPOSV, from which it will emerge
* as the inverse of S_k. 
	  


       do i=1,p
          setint(i)=i
       end do


	  i=0
	  do j=1,p
	    if(setk(j)) then
		  i=i+1
       	          setint(i)=j
	    end if
	  end do

	  do i=1,k-1
	    do j=i+1,k
                  skinput(i,j)=s(setint(i),setint(j))
		  skinput(j,i)=skinput(i,j)
	          sk(i,j)=0
                  sk(j,i)=0
            end do
            skinput(i,i)=s(setint(i),setint(i))
            sk(i,i)=1
	  end do
	  skinput(k,k)=s(setint(k),setint(k))
          sk(k,k)=1

*          dm=300
           laux='L'
           info=0
           call dposv(laux,k,k,skinput,k,sk,k,info)

* computes the trace of (inv(S_k) x (S**2)_k)**2

          dobjrv=0
          do i=1,k-1
            do j=i+1,k
              soma1=0
              soma2=0
	      do l=1,k
	        soma1=soma1+sk(i,l)*sq(setint(l),setint(j))
	        soma2=soma2+sk(j,l)*sq(setint(l),setint(i))
	      end do
	      dobjrv=dobjrv+soma1*soma2
            end do
          end do
          dobjrv=dobjrv*2
          soma=0
          do i=1,k
              soma=0
              do l=1,k
	 soma=soma+sk(i,l)*sq(setint(l),setint(i))
              end do
              dobjrv=dobjrv+soma**2
           end do
       return
       end   
********************************************************************


********************************************************************
      function dobjgcd(nqsi,qsi,valp,vecp,k,setk,p,s,fica,poriginal)
**********************************************************************
* This function computes the "variable part" of the GCD criterion. By "variable
* part" is meant that which changes from one given k-subset to the other 
* (the final value of GCD is only computed in the main subroutine, for 
* the selected solutions, in order to spare some computations in this
* subroutine, which is called extensively).  Specifically, this subroutine
* calculates the trace of the matrix (Sqsi_k) x inv(S_k), where S_k 
* is the submatrix of the covariance (or correlation) matrix for the
* variables in subset setk, "inv" stands for matrix inverse, and Sqsi_k
* 
* WARNING : This function calls the LAPACK routine DPOSV.
*
*
*  OUTPUT:
*    sum_{m in qsi} (valp_m vecp_m^t inv(S_k) vecp_m), where S_k is
*    the submatrix of the covariance/correlation matrix relevant for the 
*    k variables in subset setk, "inv" stands for matrix inverse, valp_m
*    stands for the eigenvalue of the full covariance/correlation 
*    matrix S, associated with the m-th rank (as defined in the set qsi), and
*    vecp_m is the corresponding SUB-eigenvector, defined as the subvector
*    of the eigenvector associated with valp_m, whose elements correspond to
*    positions of the k variables in setk.
*
*  INPUT: 
*
*  nqsi   - integer variable, indicating the dimension of the Principal 
*           Subspace with which the subspace spanned by the k-subset is
*            to be compared.
*  qsi    - integer vector, giving the ranks of the Principal Components 
*           which span the Principal Subspace described above (see nqsi).
*  valp   - double precision vector giving the eigenvalues of s, in 
*            decreasing order.
*  vecp   - double precision array (pxp) giving the eigenvectors of s, in
*            the same column-order as the values of valp.
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables. 
*  fica   - integer vector giving the original variable numbers of the 
*            admissible variables (those that are not forcefully excluded
*            from the subset by the user).
* poriginal -

      external dposv

      integer p,poriginal,setint(p),qsi(p),info
* integer dm
      logical setk(poriginal)
      character*1 laux
      integer fica(0:p)
      double precision s(poriginal,poriginal),skinput(k,k),sk(k,k)
      double precision valp(poriginal),vecp(poriginal,poriginal)
      double precision dobjgcd,aux0,aux


* initilizations and call to LAPACK routine DPSOV
* matrix skinput will store the submatrix S_k. Matrix s will initially
* be identity, fed to the LAPACK routine DPOSV, from which it will emerge
* as the inverse of S_k.
 	  

       do i=1,p
          setint(i)=i
       end do

	  i=0
	  do j=1,p
	    if(setk(j)) then
		  i=i+1
   	      setint(i)=j
	    end if
	  end do

	  do i=1,k-1
	    do j=i+1,k
                  skinput(i,j)=s(setint(i),setint(j))
		  skinput(j,i)=skinput(i,j)
	          sk(i,j)=0
                  sk(j,i)=0
            end do
            skinput(i,i)=s(setint(i),setint(i))
            sk(i,i)=1
	  end do
	  skinput(k,k)=s(setint(k),setint(k))
          sk(k,k)=1

*          dm=300
           laux='L'
           info=0
           call dposv(laux,k,k,skinput,k,sk,k,info)

* calculates sum_{m in qsi} (valp_m vecp_m^t inv(S_k) vecp_m)

          dobjgcd=0
	  do m=1,nqsi
	    aux=0
	    do i=1,k
		  aux0=0
		  do l=1,k
		    aux0=aux0+sk(i,l)*vecp(fica(setint(l)),qsi(m))
		  end do
                  aux0=aux0*vecp(fica(setint(i)),qsi(m))
    	          aux=aux+aux0	  
	    end do
            dobjgcd=dobjgcd+aux*valp(qsi(m))
	  end do
	  return
	  end   
************************************************************************
       function dobjtau2(k,setk,p,s,poriginal,h,rh)
**********************************************************************
* This function computes the tau_2 criterion.
* WARNING : This function calls the LAPACK routine DSYGV.
*
*  OUTPUT:
*
*   The tau^2 index of "effect magnitude".  
*   This criterion  is equivalent to the minimization of 
*   Wilks lambda statistic. tau_2 = 1 - lamlba^(1/l) Where l=min(rh,k)
*   lambda = det(E)/det(T) Where E=T-H.
*   We can also obtain lambda = d1*d2*...*dl, where d1, d2, ... dl are the 
*   eigen values of ET^-1
*  INPUT: 
*
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables.
* poriginal - 
*   h     - double precision matrix, 
*   rh    - integer, h matrix rank
*

       

       integer p,poriginal,setint(p)
         integer rh,k,lwork,itype
        integer rhaux
       logical setk(poriginal)
       character*1 laux,jobz
       double precision work(6*k),valp(k)
	 double precision s(poriginal,poriginal),h(poriginal,poriginal)
       double precision skinput(k,k),dobjtau2,lambda
	 double precision hkinput(k,k),ekinput(k,k)
        
	  
		
		external dsygv
* initializations and call to LAPACK routine dsygv
* matrix skinput will store the submatrix S_k. 
* matrix ekinput will store the submatrix E_K = S_k - H_k. 


       do i=1,p
          setint(i)=i
       end do



        i=0
       do j=1,p
         if(setk(j)) then
           i=i+1
           setint(i)=j
         end if
       end do

       do i=1,k-1
           do j=i+1,k
              skinput(i,j)=s(setint(i),setint(j))
		    hkinput(i,j)=h(setint(i),setint(j))
              skinput(j,i)=skinput(i,j)
              hkinput(j,i)=hkinput(i,j)
            end do
           skinput(i,i)=s(setint(i),setint(i))
		 hkinput(i,i)=h(setint(i),setint(i))
        end do
          skinput(k,k)=s(setint(k),setint(k))
          hkinput(k,k)=h(setint(k),setint(k))
        
	  do i=1,k-1
           do j=i+1,k
              ekinput(i,j)=skinput(i,j)-hkinput(i,j)
	        ekinput(j,i)=ekinput(i,j)
            end do
           ekinput(i,i)=skinput(i,i)-hkinput(i,i)
	   end do
	   ekinput(k,k)=skinput(k,k)-hkinput(k,k)

*          dm=300
          itype=1
		jobz='N'
		laux='L'
           info=0
	   lwork=6*k
          call dsygv(itype,jobz,laux,k,ekinput,k,skinput,k,valp,work,
     +	lwork,info)
	  
	  
	   lambda=1
	   do i=1,k
	      lambda=lambda*valp(i)
         end do
	  
	   rhaux=rh
	 if (rh.GT.(k)) then 
	    rhaux=k
	 end if 

	    dobjtau2= 1 - lambda ** (1./rhaux)


	  return
	 end   
************************************************************************
       function dobjxi2(k,setk,p,s,poriginal,h,rh)
**********************************************************************
* This function computes the xi_2 criterion.
* WARNING : This function calls the LAPACK routine DSYGV.
*
*  OUTPUT:
*
*   The criterion xi_2 of "effect magnitude".  
*   This criterion  is equivalent to the maximization of Bartlett-Pillai 
*   trace statistic U = tr(HT^-1),  xi_2 = U/l, Where l=min(rh,k)
*  
*
* INPUT: 
*
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables.
* poriginal - 
*   h     - double precision matrix
*   rh    - integer, h matrix rank
*

       

       integer p,poriginal,setint(p)
         integer rh,k,lwork,itype
       integer rhaux
       logical setk(poriginal)
       character*1 laux,jobz
       double precision work(6*k),valp(k)
	 double precision s(poriginal,poriginal),h(poriginal,poriginal)
       double precision skinput(k,k),dobjxi2,hkinput(k,k)
        
       
		
		external dsygv
* initializations and call to LAPACK routine dsygv
* matrix skinput will store the submatrix S_k. 
* matrix hkinput will store the submatrix H_k. 


       do i=1,p
          setint(i)=i
       end do

        i=0
       do j=1,p
         if(setk(j)) then
           i=i+1
           setint(i)=j
         end if
       end do

       do i=1,k-1
           do j=i+1,k
              skinput(i,j)=s(setint(i),setint(j))
		    hkinput(i,j)=h(setint(i),setint(j))
              skinput(j,i)=skinput(i,j)
              hkinput(j,i)=hkinput(i,j)
            end do
           skinput(i,i)=s(setint(i),setint(i))
		 hkinput(i,i)=h(setint(i),setint(i))
        end do
          skinput(k,k)=s(setint(k),setint(k))
          hkinput(k,k)=h(setint(k),setint(k))

*          dm=300
          itype=1
		jobz='N'
		laux='L'
           info=0
	   lwork=6*k
          call dsygv(itype,jobz,laux,k,hkinput,k,skinput,k,valp,work,
     +	lwork,info)
         dobjxi2=0
	   do i=1,k
	      dobjxi2=dobjxi2+valp(i)
         end do
	   
	   rhaux=rh
	 if (rh.GT.(k)) then 
	    rhaux=k
	 end if 
	   
	   dobjxi2=dobjxi2/rhaux  



	  return
	 end   

************************************************************************
       function dobjzeta2(k,setk,p,s,poriginal,h,rh)
**********************************************************************
* This function computes the zeta_2 criterion.
* WARNING : This function calls the LAPACK routine DSYGV.
*
*  OUTPUT:
*
*   The criterion zeta_2 of "effect magnitude".  
*   This criterion  is equivalent to the maximization of Hotteling-Lawley
*   trace statistic V = tr(HE^-1),  zeta_2 = V/(V+l), Where l=min(rh,k)
*
*  INPUT: 
*
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables.
* poriginal - 
*   h     - double precision matrix.
*   rh    - integer, h matrix rank
*

       

       integer p,poriginal,setint(p)
         integer rh,k,lwork,itype
       integer rhaux
       logical setk(poriginal)
       character*1 laux,jobz
       double precision work(6*k),valp(k)
	 double precision s(poriginal,poriginal),h(poriginal,poriginal)
       double precision skinput(k,k),dobjzeta2,hkinput(k,k)
        
       
		
		external dsygv
* initializations and call to LAPACK routine dsygv
* matrix skinput will store the submatrix S_k. 
* matrix hkinput will store the submatrix H_k.


       do i=1,p
          setint(i)=i
       end do

        i=0
       do j=1,p
         if(setk(j)) then
           i=i+1
           setint(i)=j
         end if
       end do

       do i=1,k-1
           do j=i+1,k
              skinput(i,j)=s(setint(i),setint(j))-h(setint(i),setint(j))
		    hkinput(i,j)=h(setint(i),setint(j))
              skinput(j,i)=skinput(i,j)
              hkinput(j,i)=hkinput(i,j)
            end do
           skinput(i,i)=s(setint(i),setint(i))-h(setint(i),setint(i))
		 hkinput(i,i)=h(setint(i),setint(i))
        end do
          skinput(k,k)=s(setint(k),setint(k))-h(setint(k),setint(k))
          hkinput(k,k)=h(setint(k),setint(k))
         
              
*          dm=300
          itype=1
		jobz='N'
		laux='L'
           info=0
	   lwork=6*k
          call dsygv(itype,jobz,laux,k,hkinput,k,skinput,k,valp,work,
     +	lwork,info)
         dobjzeta2=0
	   do i=1,k
	      dobjzeta2=dobjzeta2+valp(i)
         end do
	   
	   rhaux=rh
	 if (rh.GT.(k)) then 
	    rhaux=k
	 end if 
	   
	   dobjzeta2=dobjzeta2/(dobjzeta2+rhaux)  



	  return
	 end   

************************************************************************
       function dobjccr12(k,setk,p,s,poriginal,h,rh)
**********************************************************************
* This function computes the ccr1_2 criterion.
* WARNING : This function calls the LAPACK routine DSYGV.
*
*  OUTPUT:
*
*   The criterion ccr1_2 
*   This criterion  is equivalent to the maximization of Roy first
*   root (lambda1), first eigen value of HE^-1, where E=T-H.
*   ccr1_2 = lambda1/(1+lambda1). 
*  INPUT: 
*
*   k     - integer, cardinality of subset setk. 
* setk    - logical vector, indicating a given subset of variables. 
*           setk(j)=.true. iff variable j belongs to the subset. Changed
*           on output.
*   p     - integer variable, number of original variables.
*   s     - double precision matrix of covariances of all p variables.
*   h     - double precision matrix, 
*   rh    - integer
* poriginal - 


       

       integer p,poriginal,setint(p)
         integer rh,k,lwork,itype
*      integer dm
       logical setk(poriginal)
       character*1 laux,jobz
       double precision work(6*k),valp(k)
	 double precision s(poriginal,poriginal),h(poriginal,poriginal)
       double precision skinput(k,k),dobjccr12,hkinput(k,k)
        
      	external dsygv
* initializations and call to LAPACK routine dsygv
* matrix skinput will store the submatrix S_k. 
* matrix hkinput will store the submatrix H_k.


       do i=1,p
          setint(i)=i
       end do

        i=0
       do j=1,p
         if(setk(j)) then
           i=i+1
           setint(i)=j
         end if
       end do

       do i=1,k-1
           do j=i+1,k
              skinput(i,j)=s(setint(i),setint(j))-h(setint(i),setint(j))
		    hkinput(i,j)=h(setint(i),setint(j))
              skinput(j,i)=skinput(i,j)
              hkinput(j,i)=hkinput(i,j)
            end do
           skinput(i,i)=s(setint(i),setint(i))-h(setint(i),setint(i))
		 hkinput(i,i)=h(setint(i),setint(i))
        end do
          skinput(k,k)=s(setint(k),setint(k))-h(setint(k),setint(k))
          hkinput(k,k)=h(setint(k),setint(k))
         
              
*          dm=300
          itype=1
		jobz='N'
		laux='L'
           info=0
	   lwork=6*k
          call dsygv(itype,jobz,laux,k,hkinput,k,skinput,k,valp,work,
     +	lwork,info)
      
	 
	      dobjccr12=valp(k)
        
	   dobjccr12=dobjccr12/(1+dobjccr12)  



	  return
	 end   
