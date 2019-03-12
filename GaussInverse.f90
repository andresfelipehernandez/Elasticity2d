!      Subroutine GINV(n,a,x,Isym,Iwlpvt,l,u,det,Istop)
Subroutine GINV(n,a,x,Isym,Iwlpvt,det,Istop)


!  Computation of the inverse of a matrix
!  by Gauss elimination,
!  with the option of row pivoting.
! 
!
!  This subroutine returns the inverse,
!  upperand lower triangular factors,
!  flag for completion, and
!  the determinant.
!
!  SYMBOLS:
!  -------
!
!  a ...... square matrix
!  n ...... size (rows/columns) of matrix a
!  x ...... inverse of a
!
!  Isym ... flag for symmetry of matrix a (1 = symmetric)
!  Iwlpvt.. 0 for no pivoting, 1 for pivoting
!
!  eps..... tolerance to identify a singular matrix
!  tol..... tolerance for the residuals
!
!  l ...... lower triangular matrix
!  u ...... upper triangular matrix
!  det .... determinant (det(a)=det(l)*det(u))
!  Istop... flag: Istop = 1 if something is wrong
!
!  pivot .. absolute value of pivot candidates
!  ipv..... location of pivotal element
!  icount . counter for number of row interchanges
!
!---------------------------------------------------------------
integer n,na,n1,ntot,ma,m1,Istop,Icount
integer Iwlpvt,ipv,Isym
real*8 a(n,n),c(n,n+n)!,u(n,n)
real*8 l(n,n),sum,pivot,det
real*8 x(n,n),t1
real*8 eps,tol

eps=1d-8 ; tol=1d-8

!---
! initialize
!---

Istop  = 0
Icount = 0

!---
! prepare
!---
na=Isym
na = n-1
n1 = n+1
ntot = n+n
	
!---
! Initialize l and and the extended matrix c
!---

Do i=1,n
  Do j=1,n
    c(i,j) = a(i,j)
  EndDo
  Do j=1,n
    c(i,n+j) = 0.0
  EndDo
  c(i,n+i) = 1.0
EndDo

!---
!  begin row reductions
!---
	
   Do m = 1,na            ! outer loop for working row

      ma = m-1
      m1 = m+1

      If(Iwlpvt.ne.1) Go to 97    ! skip pivoting

!-----------------------------
!  Pivoting
!  begin by searching column i
!  for largest element
!-----------------------------
	   
      ipv  = m
      pivot = dabs(c(m,m))
	 
      Do j = m1,n
        If(dabs(c(j,m)).gt.pivot) then
	   ipv   = j
         pivot = dabs(c(j,m))
        EndIf
    EndDo

      If(pivot.lt.eps) then
        ! write (6,*) 
        ! write (6,*) " Trouble at station 1 of Gauss elimination"
        ! write (6,*) 
        Istop = 1
        Return
      EndIf

!---
!  switch the working row with the row containing the 
!  pivot element (also switch rows in l)
!---
	 
      If(ipv.ne.m) then

        Do j = m,ntot
		t1    = c(m,j)
          c(m,j)  = c(ipv,j)
          c(ipv,j)= t1
      EndDo

        Do j=1,ma
         t1     = l(m,j)
         l(m,j)   = l(ipv,j)
         l(ipv,j) = t1
        EndDo

	  Icount=Icount+1 
	
      EndIf

97    Continue

!---------------------------------------
! reduce column i beneath element !(m,m)
!---------------------------------------
	 
      Do  i = m1,n

        If(Isym.eq.1) then

          l(i,m) = c(m,i)/c(m,m)
          Do j = i,ntot
           c(i,j)=c(i,j)-l(i,m)*c(m,j)
          EndDo

        Else

          l(i,m) = c(i,m)/c(m,m)
          c(i,m) = 0.0
	  Do j=m1,ntot
	   c(i,j)=c(i,j)-l(i,m)*c(m,j)
      EndDo

	  EndIf

    enddo

	
!---
!  fill in the leading zeros in the matrix c
!  (not necessary)
!---
	 
      Do j=1,m
        c(m+1,j)=0.0
      EndDo

    enddo         ! end of outer loop for working row

!---
! check the last diagonal element
!  for singularity
!---
	
      If(dabs(c(n,n)).lt.eps) then
        ! write (6,*)
        ! write (6,*) " Trouble at station 2 of Gauss elimin."
        Istop = 1
        Return
      EndIf

!---
! complete the matrix l
!---

    Do  i=1,n
        l(i,i)=1.0
    EndDo
	
!---
! define the matrix u
!---

!      Do 71 i =1,n
!        Do 72 j =1,n
!         u(i,j) = c(i,j)
!72      EndDo
!71    EndDo

!---
! perform back-substitution to solve
! the reduced system
! using the upper triangular matrix c
!---
	
      Do  ll = 1,n

        x(n,ll) = c(n,n+ll)/c(n,n)

	  Do  i=n-1,1,-1
	  sum=c(i,n+ll)
	  Do  j=i+1,n
	   sum=sum-c(i,j)*x(j,ll)
      EndDo
	  x(i,ll) = sum/c(i,i)
      EndDo

    EndDo

!---
! compute determinant as:
! det(a)=det(l)*det(u)
!---
	
      det=1.0
	
      Do  i=1,n
        det=det*c(i,i)
    EndDo

      If(Iwlpvt.eq.1) then
        ! write (6,*)
        ! write (6,*) " Number of row interchanges : ",Icount
        ! write (6,*)
        Do  i=1,icount
          det = -det
      EndDo
      EndIf

!---
! compute residuals
!---

      Do   ll=1,n

       Do   i =1,n

        sum = 0.0
        If(i.eq.ll) sum = 1.0

        Do  j=1,n
         sum = sum - a(i,j)*x(j,ll)
      EndDo

        If(dabs(sum).gt.tol) then
          Istop = 1
          ! write (6,*) "gel_inv: failed to compute the inverse"
          !write (6,100) i,sum
        EndIf

     EndDo
    EndDo

!---
! wrap up
!---
!100   Format (1x,i4,f15.10)
!101   Format(16(16(1x,f5.3),/))

Return
End
