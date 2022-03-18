!*************************************************************************
!
! Module to compute the maxent basis functions using Newton's method
! with backtracking line search algorithm
! Programmer: Subhajit Banerjee (jitbanerjee@ucdavis.edu)
! November, 2014
! V1.0
! ++++++++
! Purpose:
! ++++++++
! 
! +++++++++++++
! Dependencies:
! +++++++++++++
!
!*************************************************************************

MODULE MAXENT

USE NRTYPE

IMPLICIT NONE
CONTAINS

!*************************************************************************

SUBROUTINE  jaynesdual(fval, dfval, ddfval, xa, x, lambda, nsdim, n, W)
! Computes the Objective function for the dual problem when the primal
! problem uses Jaynes relative entropy as the basis function
! NOTE: Primal proble is -- Maximize -(KL distance) + constaints
! Hence, Dual problem is -- Minimize (logZ)
! Interfacing Variables:
INTEGER(I4B), INTENT(IN)                                  :: nsdim, n
REAL(DP), INTENT(OUT)                                     :: fval
REAL(DP), DIMENSION(nsdim), INTENT(OUT), OPTIONAL         :: dfval
REAL(DP), DIMENSION(nsdim,nsdim), INTENT(OUT), OPTIONAL   :: ddfval
REAL(DP), DIMENSION(nsdim), INTENT(IN)                    :: x, lambda
REAL(DP), DIMENSION(n,nsdim)                              :: xa
REAL(DP), DIMENSION(n),  INTENT(IN)                       :: W

! Local variables:
REAL(DP), DIMENSION(n)                          :: Za
INTEGER(I4B)                                    :: a, b, i, j
REAL(DP)                                        :: temp, Z

fval = 0.d0
dfval = 0.d0
ddfval = 0.d0

! Objective function = log(Z)
do a = 1,n
    Za(a) =  W(a)*exp(-dot_product(lambda,(xa(a,:)-x)))
enddo

Z = sum(Za)
fval = log(Z)             ! for minimization

! Gradient of objective function w.r.t. \lambda, i.e.,
! \nabla_{\lambda} [log(Z)]
if (PRESENT(dfval)) then
    do i=1,nsdim
        do a = 1,n
            dfval(i) =  dfval(i) + Za(a)*(xa(a,i)-x(i))
        enddo
    enddo
    dfval = -dfval/Z      ! for minimization 
endif

! Hessian of objective function w.r.t. \lambda, i.e.,
! \nabla_{\lambda}*\nabla_{\lambda} [log(Z)]
if (PRESENT(ddfval)) then
    do i=1,nsdim
        do j =1,nsdim

            temp = 0.d0
            do b=1,n
                temp = temp + Za(b)*(xa(b,j)-x(j))
            enddo
            temp = temp/Z

            do a=1,n
                ddfval(i,j) = ddfval(i,j) + (Za(a)*(xa(a,i) - x(i))* &
                    (xa(a,j)- x(j))/Z - Za(a)*(xa(a,i) - x(i))*temp/Z)
            enddo
        enddo
    enddo
endif

END SUBROUTINE jaynesdual

!*************************************************************************
! Updated 10/24/2001.
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DTRM (A,N,D,INDX)
!
! Subroutine for evaluating the determinant of a matrix using 
! the partial-pivoting Gaussian elimination scheme.
! Copyright (c) Tao Pang 2001.
!
!  IMPLICIT NONE
INTEGER(I4B), INTENT (IN)                 :: N
INTEGER(I4B)                              :: I,J,MSGN
INTEGER(I4B), INTENT (OUT), DIMENSION (N) :: INDX
REAL(DP), INTENT (OUT)                    :: D
REAL(DP), INTENT (INOUT), DIMENSION (N,N) :: A
!
CALL ELGS(A,N,INDX)
!
D = 1.d0
DO I = 1, N
    D = D*A(INDX(I),I)
ENDDO
!
MSGN = 1
DO I = 1, N
    DO WHILE (I.NE.INDX(I))
        MSGN = -MSGN
        J = INDX(I)
        INDX(I) = INDX(J)
        INDX(J) = J
    ENDDO
ENDDO
D = MSGN*D
END SUBROUTINE DTRM
!
SUBROUTINE ELGS (A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
!IMPLICIT NONE
INTEGER(I4B), INTENT (IN)                 :: N
INTEGER(I4B)                              :: I,J,K,ITMP
INTEGER(I4B), INTENT (OUT), DIMENSION (N) :: INDX
REAL(DP)                                  :: C1,PII,PII1,PJ
REAL(DP), INTENT (INOUT), DIMENSION (N,N) :: A
REAL(DP), DIMENSION (N)                   :: C
!
! Initialize the index
!
DO I = 1, N
    INDX(I) = I
ENDDO
!
! Find the rescaling factors, one from each row
!
DO I = 1, N
    C1= 0.d0
    DO J = 1, N
      C1 = AMAX1(C1,ABS(A(I,J)))
    ENDDO
    C(I) = C1
ENDDO
!
! Search the pivoting (largest) element from each column
!
DO J = 1, N-1
    PII1 = 0.d0
    DO I = J, N
      PII = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PII.GT.PII1) THEN
        PII1 = PII
        K   = I
      ENDIF
    ENDDO
!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      ENDDO
    ENDDO
ENDDO
!
END SUBROUTINE ELGS

!*************************************************************************

SUBROUTINE invmatrix(A, invA, N)
! Interfacing Variables:
REAL(DP), DIMENSION(N,N), INTENT(OUT)              :: invA
REAL(DP), DIMENSION(N,N), INTENT(IN)               :: A
INTEGER(I4B), INTENT(IN)                           :: N

! Local Variables:
INTEGER(I4B)                                       :: M, LDA, INFO, LWORK
INTEGER, DIMENSION(N)                              :: IPIV
REAL(DP), ALLOCATABLE, DIMENSION(:)                :: WORK

M = N
LDA = N

! Store A in Ainv to prevent it from being overwritten by LAPACK
invA = A
! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
CALL DGETRF(M, N, invA, LDA, IPIV, INFO)

IF (INFO.EQ.0) THEN
!    PRINT '(" LU decomposition successful ")'
    CONTINUE
ENDIF

IF(INFO.LT.0)THEN
    PRINT '(" LU decomposition:  illegal value ")'
    STOP
ENDIF

IF(INFO.GT.0)THEN
    WRITE(*,35)INFO,INFO
35  FORMAT( 'LU decomposition: U(',I4,',',I4,') = 0 ')
ENDIF

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
LWORK = N*N
ALLOCATE ( WORK(LWORK) )
CALL DGETRI(N, invA, N, IPIV, WORK, LWORK, INFO)

IF (info.NE.0) THEN
    stop 'Matrix inversion failed!'
ENDIF

END SUBROUTINE invmatrix

!*************************************************************************

SUBROUTINE newton(lambda_min, max_entropy_min, lambda_iter, &
                  nsdim, n, objectivefuncname, xa, x, W, MAXITER, TOL)

USE computeprior

! Interfacing Variables:
INTEGER(I4B), INTENT(IN)                            :: MAXITER
REAL(DP)                                            :: TOL
INTEGER(I4B), INTENT(IN)                            :: n, nsdim
REAL(DP), DIMENSION(MAXITER,nsdim), INTENT(OUT)     :: lambda_iter
REAL(DP), DIMENSION(nsdim), INTENT(OUT)             :: lambda_min
REAL(DP), INTENT(OUT)                               :: max_entropy_min
CHARACTER(80), INTENT(IN)                           :: objectivefuncname
REAL(DP), DIMENSION(nsdim), INTENT(IN)              :: x
REAL(DP), DIMENSION(n,nsdim)                        :: xa
REAL(DP), DIMENSION(n), INTENT(IN)                  :: W

! Local Variables:
REAL(DP), DIMENSION(nsdim)                          :: lambda ! Initial guess
REAL(DP), DIMENSION(nsdim)                          :: lambdanew
INTEGER(I4B)                                        :: i, j, k
REAL(DP)                                            :: alphak
REAL(DP), PARAMETER                                 :: RHO = 0.5d0
REAL(DP), PARAMETER                                 :: C1 = 1.d-04
REAL(DP), DIMENSION(MAXITER,nsdim)                  :: searchdir
REAL(DP), DIMENSION(MAXITER)                        :: step
REAL(DP), DIMENSION(nsdim,nsdim)                    :: ddfval, ddfvalnew
REAL(DP), DIMENSION(nsdim,nsdim)                    :: invddfval
REAL(DP), DIMENSION(nsdim)                          :: dfval, dfvalnew, pk
REAL(DP)                                            :: fval, fvalnew

lambda = 0.0d0

select case (objectivefuncname)
    case ("jaynes")
        ! Call the corresponding function from computeprior module:
        ! Objec. func. call # 1, with derivatives:
        CALL jaynesdual(fval, dfval, ddfval, xa, x, lambda, nsdim, n, W)
    case default
        write(*,*)"Coming Soon! Exiting..."
        STOP
end select

do k = 1,MAXITER
    ! Termination condition:
    if ( sqrt(dot_product(dfval,dfval)) <= TOL ) exit

    ! Newton direction:
    CALL invmatrix(ddfval, invddfval, nsdim)
    pk = -matmul(invddfval,dfval)

    ! Steepest descent direction:
    ! pk = -dfval
    alphak = 1.d0
    
    ! Updating the design variables, \lambda:
    lambdanew = lambda + alphak*pk      ! Guess value

    ! Objec. func. call # 2, with derivatives:
    CALL jaynesdual(fvalnew, dfvalnew, ddfvalnew, xa, x, lambdanew, nsdim, n, W)

    ! Backtracking line search algorithm:
    if (fvalnew > (fval + C1*alphak*dot_product(dfval,pk))) then ! NOT Sufficient decrease
        alphak = alphak*RHO
    endif

    step(k) = alphak
    lambdanew = lambda + alphak*pk      ! Value used to compute next step

    ! Objec. func. call # 3, with derivatives:
    CALL jaynesdual(fvalnew, dfvalnew, ddfvalnew, xa, x, lambdanew, &
                    nsdim, n, W)

    lambda_iter(k,:) = lambdanew
    searchdir(k,:) = pk

    lambda = lambdanew
    dfval = dfvalnew
    ddfval = ddfvalnew
    fval = fvalnew

    if (k == MAXITER) then
        write(*,*)"ERROR!"
        write(*,*)"The algorithm did not converge in", MAXITER, "iterations."
        write(*,*)"EXITING..."
        STOP
    endif

enddo

lambda_min = lambda
max_entropy_min = fval

END SUBROUTINE newton

!*************************************************************************

SUBROUTINE shapefunc(phi,dphi,n,nsdim,lambda,weightfuncname,Wa,dWa,fval,xa,x)

! Interfacing Variables:
REAL(DP), DIMENSION(n), INTENT(OUT)         :: phi       ! max-ent basis
REAL(DP), DIMENSION(n,nsdim), INTENT(OUT)   :: dphi      ! grad(max-ent)
INTEGER(I4B), INTENT(IN)                    :: n, nsdim
REAL(DP), DIMENSION(nsdim), INTENT(IN)      :: lambda    ! Optimum lambda from Dual problem
REAL(DP), DIMENSION(nsdim), INTENT(IN)      :: x
REAL(DP), DIMENSION(n,nsdim), INTENT(IN)    :: xa
CHARACTER(80), INTENT(IN)                   :: weightfuncname
REAL(DP), DIMENSION(n), INTENT(IN)          :: Wa        ! Nodal prior weight function
REAL(DP), DIMENSION(n,nsdim), INTENT(IN)    :: dWa       ! deriv. prior weight function
REAL(DP), INTENT(IN)                        :: fval      ! Minimized [log(Z)]

! Local Variables:
INTEGER(I4B)                            :: a, b, i, j
REAL(DP)                                :: Z, temp
REAL(DP), DIMENSION(n)                  :: fa, ga, Za
REAL(DP), DIMENSION(nsdim,nsdim)        :: Jopti, Dlambdaopti
REAL(DP), DIMENSION(nsdim,nsdim)        :: invJopti
REAL(DP), DIMENSION(nsdim)              :: gbdwb
REAL(DP), DIMENSION(nsdim,nsdim)        :: Hopti, invHopti

! All are computed at the optimum solution \lambda^*
! Ref: Cell-based max-ent paper:

! Z = 
Z = exp(fval)
! phi =
do a =1,n
    ! Approach # 1:
    fa(a) = exp( dot_product(lambda, (x - xa(a,:))) )
    ga(a) = fa(a)/Z
    phi(a) = Wa(a)*ga(a)
enddo

select case (weightfuncname)
    case ("gaussian-rbf")       ! Short and sweet expression for Gaussian-rbf
        do a =1,n
            Za(a) = Wa(a)*fa(a)
        enddo

        Hopti = 0.d0
        do i=1,nsdim
            do j =1,nsdim
                do a=1,n
                    Hopti(i,j) = Hopti(i,j) + Za(a)*(xa(a,i) - x(i))* &
                                (xa(a,j)- x(j))/Z
                enddo
            enddo
        enddo
        CALL invmatrix(Hopti, invHopti, nsdim)
        do a=1,n
            dphi(a,:) = phi(a)*matmul(invHopti, (xa(a,:)-x))
        enddo
    case default
        ! J* = Checked
        Jopti = 0.d0
        do i=1,nsdim
            do j=1,nsdim
                do a=1,n
                    Jopti(i,j) = Jopti(i,j) + phi(a)*(x(i)-xa(a,i))* &
                                            (x(j)-xa(a,j))
                enddo
            enddo
        enddo

        ! D \lambda* = 
        Dlambdaopti = 0.d0
        do i=1,nsdim
            do j=1,nsdim
                do a=1,n
                    if (i .EQ. j) then
                        Dlambdaopti(i,j) = Dlambdaopti(i,j) + &
                        ga(a)*dWa(a,i)*(x(j)-xa(a,j)) + 1.d0/n
                    elseif (i .NE. j) then
                         Dlambdaopti(i,j) = Dlambdaopti(i,j) + &
                         ga(a)*dWa(a,i)*(x(j)-xa(a,j))
                    endif
                enddo
            enddo
        enddo

        CALL invmatrix(Jopti,invJopti,nsdim)
        !Dlambdaopti = -matmul(Dlambdaopti,invJopti)         ! Check!
        Dlambdaopti = -matmul(invJopti,Dlambdaopti)

        gbdwb = 0.d0
        do i=1,nsdim
            do a=1,n
                gbdwb(i) = gbdwb(i) + ga(a)*dWa(a,i)
            enddo
        enddo

        do a=1,n
    !    dphi(a,:) = ga(a)*dWa(a,:) + phi(a)* &
    !            (matmul(Dlambdaopti, (x-xa(a,:))) - gbdwb)
        dphi(a,:) = ga(a)*dWa(a,:) + phi(a)* &
            (matmul((x-xa(a,:)), Dlambdaopti) - gbdwb)

        enddo
end select

END SUBROUTINE shapefunc

!*************************************************************************

SUBROUTINE driver(n, nsdim, weightfuncname, epsilon0, & 
        objectivefuncname, Xa, X, Ra, MAXITER, eps, phi, dphi)

USE COMPUTEPRIOR
INTEGER(I4B), INTENT(IN)            :: n, nsdim
REAL(DP), INTENT(IN)                :: epsilon0, Ra
CHARACTER(80), INTENT(IN)           :: weightfuncname, objectivefuncname
REAL(DP), INTENT(IN)                :: Xa(n,nsdim)
REAL(DP), INTENT(IN)                :: X(nsdim)
INTEGER(I4B), INTENT(IN)            :: MAXITER
REAL(DP), INTENT(IN)                :: eps
REAL(DP), INTENT(OUT)               :: phi(n)
REAL(DP), INTENT(OUT)               :: dphi(n,nsdim)

! Local Variables:
INTEGER(I4B)                        :: i, j
! For Objective Function & Weight Function:
REAL(DP), DIMENSION(n)              :: Wa
REAL(DP), DIMENSION(n,nsdim)        :: dWa
! For Newton's Scheme:
REAL(DP), DIMENSION(MAXITER,nsdim)  :: lambda_iter
REAL(DP), DIMENSION(nsdim)          :: lambda_min
REAL(DP)                            :: max_entropy_min

select case (weightfuncname)
    case ("gaussian-rbf")
        ! Call the corresponding function from computeprior module:
        CALL gaussianrbf(Wa, dWa, x, xa, n, nsdim, Ra, epsilon0)
    case default
        write(*,*)"Coming Soon! Exiting..."
end select

CALL newton(lambda_min, max_entropy_min, lambda_iter, &
                  nsdim, n, objectivefuncname, xa, x, Wa, MAXITER, eps)


CALL shapefunc(phi, dphi, n, nsdim, lambda_min, weightfuncname, Wa, dWa, &
    max_entropy_min, xa, x)

END SUBROUTINE driver

END MODULE maxent
