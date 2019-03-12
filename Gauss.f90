      SUBROUTINE SLNPD(A,XT,B,N)
!
!  PROGRAM 6 
! 
! SOLUTION OF LINEAR SYSTEMS OF EQUATIONS 
! BY THE GAUSS ELIMINATION METHOD PROVIDING 
! FOR INTERCHANGING ROWS WHEN ENCOUNTERING A
! ZERO DIAGONAL COEFICIENT
! 
! A : SYSTEM MATRIX 
! B : LINEAR SYSTEM SOLUTION. 
! X:  INDEPENDENT TERMS
! 
! N : ACTUAL NUMBER OF UNKNOWNS 
! NX: ROW AND COLUMN DIMENSION OF A 
! 
        INTEGER N,N1,K1,K
      REAL*8 B(N),A(N,N),XT(N,1),TOL,C

        DO 91 I=1,N
         B(I)=XT(I,1)
91      CONTINUE

!
      TOL=1.E-8
!
      N1=N-1
      DO 100 K=1,N1 
      K1=K+1
      C=A(K,K)
      IF(DABS(C)-TOL)1,1,3
1       DO 7 J=K1,N 
! 
! TRY TO INTERCHANGE ROWS TO GET NON ZERO DIAGONAL COEFFICIENT
! 
      IF(DABS((A(J,K)))-TOL)7,7,5 
5       DO 6 L=K,N
      C=A(K,L)
      A(K,L)=A(J,L) 
6     A(J,L)=C
      C=B(K)
      B(K)=B(J) 
      B(J)=C
      C=A(K,K)
      GO TO 3 
    7 CONTINUE
      GO TO 8
! 
! DIVIDE ROW BY DIAGONAL COEFFICIENT
! 
3     C=A(K,K)
      DO 4 J=K1,N 
4     A(K,J)=A(K,J)/C
      B(K)=B(K)/C
! 
! ELIMINATE UNKNOWN X(K) FROM ROW I 
! 
      DO 10 I=K1,N
      C=A(I,K)
      DO 9 J=K1,N 
9     A(I,J)=A(I,J)-C*A(K,J)
10    B(I)=B(I)-C*B(K)
100   CONTINUE
! 
! COMPUTE LAST UNKNOWN
! 
      IF(DABS((A(N,N)))-TOL)8,8,101 
101   B(N)=B(N)/A(N,N)
! 
! APPLY BACKSUBSTITUTION PROCESS TO COMPUTE REMAINING UNKNOWNS
! 
      DO 200 L=1,N1 
      K=N-L 
      K1=K+1
      DO 200 J=K1,N 
200   B(K)=B(K)-A(K,J)*B(J) 
! 
! COMPUTE VALUE OF DETERMINANT
! 
      D=1.
      DO 250 I=1,N
       D=D*A(I,I)
250     CONTINUE

      GO TO 300
8     PRINT*, "**** SINGULARITY IN ROW", K
      PAUSE
      D=0.


300   RETURN
      END 
!-----------------------------------------------------------------------  

