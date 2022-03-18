!*************************************************************************
!                MAIN PROGRAM TO TEST MAXENT MODULE                      !
!*************************************************************************
!                       Commmand to compile:                             !
!    gfortran -g -fcheck=all -Wall -o minxent minxent2dpoisson.f90       !

PROGRAM MAIN

! f90 main program: reads input, processes, writes output
USE NRTYPE
USE MAXENT
USE COMPUTEPRIOR

IMPLICIT NONE
REAL(DP), PARAMETER     :: epsilon0 = 1.d-06
INTEGER(I4B)            :: n, nsdim  ! number of nodes and space-dimension
REAL(DP), ALLOCATABLE   :: x(:)      ! point (mostly Gauss) in nsdim-dimensions
REAL(DP), ALLOCATABLE   :: xa(:,:)   ! nodal coordinates (matrix)

CHARACTER(80)           :: weightfuncname = "gaussian-rbf" ! ``prior'' in the weight function
CHARACTER(80)           :: objectivefuncname = "jaynes"
REAL(DP), ALLOCATABLE   :: phi(:)      ! basis functions
REAL(DP), ALLOCATABLE   :: dphi(:,:)   ! first derivatives of basis functions

CHARACTER(80)         :: inputfile       ! Sample input file
CHARACTER(80)         :: outputfile      ! Sample output file
INTEGER(I4B)          :: i, j, a, ok
INTEGER(I4B)          :: maxit
REAL(DP)              :: eps, h, Ra

inputfile = "sample.in"
outputfile = "sample.out"

! open sample input file 
OPEN(UNIT=21,file = inputfile,STATUS='OLD',ACTION='READ')   ! read from sample file
OPEN(UNIT=22,file = outputfile,STATUS='NEW')                ! write to sample file

read(21,*,END=100)                    ! header: NSD and NODES
read(21,*,END=100)nsdim,n             ! dimension and # of nodes
allocate ( xa(n,nsdim) )              ! Nodal Coordinates
allocate ( x(nsdim) )                 ! Point of computation within conv hull
allocate ( phi(n) )                   ! Shape functions                   
allocate ( dphi(n,nsdim) )            ! Derivatives of Shape functions

read(21,*,END=100)              ! header: POINT
read(21,*,END=100)x
read(21,*,END=100)              ! header: COORDS
do i = 1,n
  read(21,*,END=100)xa(i,:)
enddo

! h comes from nodal spacing (typical/max nodal spacing):
!h = abs(xa(2,1) - xa(1,1)) ! Valid for uniform mesh
h = 1.d0
Ra = 2.d0*h

read(21,*,END=100)              ! header: MAXITER AND TOL
read(21,*,END=100)maxit,eps

read(21,*,END=100)              ! header: PRIOR AND OBJECTIVE
read(21,*)weightfuncname, objectivefuncname

100 CLOSE(UNIT=21)

CALL driver(n, nsdim, weightfuncname, epsilon0, objectivefuncname, &
                xa, x, Ra, maxit, eps, phi, dphi)

write(*,*) "Now writing output to an external data file..."
write(*,*) "..."
write(*,*) " "
write(22,*)"Spatial dimension = ", nsdim
write(22,*)"Total number of nodes = ", n
write(22,*)" "
! Nodal Coordinates:
write(22,*)"Nodal Coordinates: "
do a = 1,n
  write(22,*)"Node = ",a," Coordinates = ", ( Xa(a,i), i = 1,nsdim )
enddo
write(22,*)" "
! Problem Setup details:
write(22,*)"Maximum number of iterations = ", maxit
write(22,*)"Convergence tolerance = ", eps
write(22,*)"Nodal Spacing = ", h
write(22,*)"Support radius = ", Ra
write(22,*)"Prior weight function: ", weightfuncname
write(22,*)"Objective Functional: ", objectivefuncname
write(22,*)" "
! Solution Details:
write(22,*)"Value of Maxent Shape Functions at Point: " 
do i=1,nsdim
    write(22,*) X(i)
enddo
write(22,*)"is: "
write(22,*) " "
! PHI
do a=1,n
    write(22,*) phi(a)
enddo
write(22,*)" "
write(22,*)"The sum of the maxent shape functions = ", sum(phi)
write(22,*)"Value of the derivatives of the shape functions: "
write(22,*)
! grad(PHI)
do a = 1,n
    write(22,*) (dphi(a,i), i = 1,nsdim )
enddo

CLOSE(UNIT=22)
END PROGRAM main
