!*************************************************************************
!
! Module to Compute the prior weight-functions as used in max-ent basis
! Programmer: Subhajit Banerjee (jitbanerjee@ucdavis.edu)
! November, 2015
! V1.9
! ++++++++
! Purpose:
! ++++++++
! 
! +++++++++++++
! Dependencies:
! +++++++++++++
!
!*************************************************************************

MODULE COMPUTEPRIOR
USE NRTYPE

IMPLICIT NONE

CONTAINS

!*************************************************************************
SUBROUTINE gaussianrbf(Wa, dWa, x, xa, n, nsdim, Ra, epsilon0)
! Interfacing Variables:
IMPLICIT NONE
REAL(DP), INTENT(IN)                        :: Ra, epsilon0
INTEGER(I4B), INTENT(IN)                    :: n, nsdim
REAL(DP), DIMENSION(n,nsdim),INTENT(IN)     :: xa
REAL(DP), DIMENSION(nsdim), INTENT(IN)      :: x
REAL(DP), DIMENSION(n), INTENT(OUT)         :: Wa
REAL(DP), DIMENSION(n,nsdim), INTENT(OUT)   :: dWa
! Local variable
REAL(DP)                                    :: beta
REAL(DP)                                    :: L2norm
INTEGER(I4B)                                :: a, i
REAL(DP)                                    :: macheps

macheps = epsilon(0d0)

beta = -log(epsilon0)/Ra**2

do a=1,n
    L2norm = sqrt(dot_product((x-xa(a,:)),(x-xa(a,:))))
    if (L2norm >= (Ra-100.d0*macheps)) then
        Wa(a) = 0.d0
    else
        Wa(a) = exp(-beta*L2norm**2)
    endif
enddo

do a=1,n
    do i=1,nsdim
        L2norm = sqrt(dot_product((x-xa(a,:)),(x-xa(a,:))))
        if (L2norm > (Ra-100.d0*macheps)) then
            dWa(a,i) = 0.d0
        else
            dWa(a,i) = -2.d0*beta*exp(-beta*L2norm**2)*(x(i) - xa(a,i))
        endif
    enddo
enddo

END SUBROUTINE gaussianrbf

!*************************************************************************

END MODULE computeprior
