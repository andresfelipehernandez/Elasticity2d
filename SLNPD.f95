SUBROUTINE SLNPD

!==================================================================================================================
!Subrutina encargada de solucionar el sistema AX=B generada en la subrutina assembly
!==================================================================================================================

use vars

integer N1,K1,L1
real C,CD,CDD

!allocate (A(n,n)),(B(n))


!       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!
! 1 SOLUTION OF THE LINEAR SYSTEM OF EQUATIONS
!  BY THE GAUSS ELIMINATION METHOD PROVIDING
!  FOR INTERCHANGING ROWS WHEN ENCOUNTERING A
!  ZERO DIAGONAL COEFICIENT
!
!  A : SYSTEM MATRIX
! DFI: ORIGINALLY IT CONTAINS THE INDEPENDENT
!      COEFFICIENTS, AFTER SOLUTION IT CONTAIN
!      THE VALUES OF THE SYSTEM UNKNOWNS
!
! KC : ACTUAL NUMBER OF UNKNOWNS

      N1=nn-1
      DO 100 K=1,N1
      K1=K+1
      C=ma(K,K)
      IF(ABS(C)-0.000000000001)1,1,3
    1 DO 7 J=K1,nn

!C  TRY TO INTERCHANGE ROWS TO GET NON-ZERO DIAGONAL COEFFICIENT

      IF(ABS(ma(J,K))-0.000000000001)7,7,5
    5 DO 6 L1=K,nn
      C=ma(K,L1)
      ma(K,L1)=ma(J,L1)
    6 ma(J,L1)=C
      C=b(K,1)
      b(K,1)=b(J,1)
      b(J,1)=C
      C=ma(K,K)
      GO TO 3
    7 CONTINUE
  778 WRITE(IMP,2) K
    2 FORMAT(' **** SINGULARITY IN ROW',I5)
      GO TO 300

!C  DIVIDE ROW BY DIAGONAL COEFFICIENT

    3 C=ma(K,K)
      DO 4 J=K1,nn
    4 ma(K,J)=ma(K,J)/C
      b(K,1)=b(K,1)/C

!C  ELIMINATE UNKNOWN X(K) FROM ROW I

      DO 10 I=K1,nn,1
      CD=ma(I,K)
      DO 9 J=K1,nn,1
      CDD=CD
      ma(I,J)=ma(I,J)-CDD*ma(K,J)
    9 CONTINUE
      b(I,1)=b(I,1)-CD*b(K,1)
   10 CONTINUE
  100 CONTINUE

!C  COMPUTE LAST UNKNOWN

      IF(ABS(ma(nn,nn))-0.000000000001) 777,777,101
  777 WRITE(IMP,22) K
   22 FORMAT(' **** SINGULARITY IN ROW',I5)
      GO TO 300
  101 b(nn,1)=b(nn,1)/ma(nn,nn)

!C  APPLY BACKSUBSTITUTION PROCESS TO COMPUTE REMAINING UNKNOWNS

      DO 200 L1=1,N1
      K=nn-L1
      K1=K+1
      DO 200 J=K1,nn
  200 b(K,1)=b(K,1)-ma(K,J)*b(J,1)
      open(30,file='Solucionsistema.txt',status='unknown')
      write(30,*)'La solucion del sistema es'
      write(30,*)B
      close(30)
  300 RETURN
      END
