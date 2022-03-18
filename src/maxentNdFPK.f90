!*************************************************************************
!                                                                        !
!                           PROGRAM MAXENT2DFPK                          !
!                           ===================                          !
!                                                                        !
!                                                                        !      
! Purpose                                                                !
! =======                                                                !
! FORTRAN 90 Program to solve two-dimensional Fokker-Planck-Kolmogorov   !
! equation using Minimum Cross Entropy (MinXEnt) basis function          !
! Written by    : Subhajit Banerjee, UC Davis                            !
! Date          : January 2015                                           !
! Dependencies  : Fortran 90 Modules --                                  !
!                1. computeprior.f90                                     !
!                2. maxent.f90                                           !
!                                                                        !
!                                                                        !
!*************************************************************************
      
PROGRAM maxent2dFPK

USE COMPUTEPRIOR        ! use types and functions defined in this module 
USE MAXENT              ! use types and functions defined in this module
USE NRTYPE

IMPLICIT NONE           ! require explicit type declaration/checking

! Interface with called subroutines:
INTERFACE
    SUBROUTINE gaussq_triangle(nint,gauss)
        USE NRTYPE
        INTEGER(I4B), INTENT(IN)    :: nint
        REAL(DP), INTENT(OUT)       :: gauss(3,nint)
    END SUBROUTINE
END INTERFACE

INTERFACE
    SUBROUTINE  egauss(gs,nsdim,nodcoord,conn,gauss,ntriangle,numnod,numq2, quado)
        USE NRTYPE
        integer(i4b), intent(in)                         :: ntriangle, numnod
        integer(i4b), intent(in)                         :: numq2, quado, nsdim
        integer(i4b), dimension(3,ntriangle), intent(in) :: conn
        real(dp), dimension(3,quado), intent(in)         :: gauss         
        real(dp), dimension(nsdim,numnod), intent(in)    :: nodcoord      
        real(dp), dimension(4,numq2), intent(out)        :: gs            
    END SUBROUTINE
END INTERFACE

INTERFACE
    SUBROUTINE domain(v,gpos,nodcoord,Ra,nnodes)
        USE NRTYPE
        real(DP), dimension(2),  intent(in)                  :: gpos
        integer(I4B), intent(in)                             :: nnodes
        real(DP), intent(in)                                 :: Ra
        real(DP), dimension(2,nnodes), intent(in)            :: nodcoord
        integer(I4B), allocatable, intent(out), dimension(:) :: v
    END SUBROUTINE
END INTERFACE

INTERFACE
    SUBROUTINE getcoordinates(xyz,v,numnghbr,nnodes,nodcoord)
        USE NRTYPE
        integer(I4B), dimension(numnghbr),  intent(in) :: v
        integer(I4B), intent(in)                       :: nnodes, numnghbr
        real(DP), dimension(2,nnodes), intent(in)      :: nodcoord
        real(DP), dimension(numnghbr,2), intent(out)   :: xyz
    END SUBROUTINE
END INTERFACE

INTERFACE
    SUBROUTINE supportplot(gs,nodcoord,Ra,numq2,nnodes,filename)
        USE NRTYPE
        integer(I4B), intent(in)                     :: numq2, nnodes
        real(DP), dimension(4,numq2), intent(out)    :: gs
        real(DP), dimension(2,nnodes), intent(in)    :: nodcoord
        real(DP), intent(in)                         :: Ra
        character(80),  intent(in)                   :: filename
    END SUBROUTINE
END INTERFACE

INTERFACE
    SUBROUTINE gridcircle(allnodes,nodloc,Ra,ndivtheta,ndivr)
        USE NRTYPE
        integer(I4B), intent(in)                                       :: ndivtheta, ndivr
        real(DP), intent(in)                                           :: Ra
        real(DP), dimension(2), intent(in)                             :: nodloc
        real(DP), dimension(2,(1 + (ndivtheta+1)*ndivr)), intent(out)  :: allnodes
    END SUBROUTINE
END INTERFACE
 
INTERFACE
    SUBROUTINE initialcond(p0, nodcoord, h, A, sigma, x0, nsdim, numnod)
        USE NRTYPE
        integer(i4b), intent(in)                        :: nsdim, numnod
        real(dp), intent(in), dimension(nsdim,numnod)   :: nodcoord
        real(dp), intent(in), dimension(nsdim)          :: x0
        real(dp), intent(in)                            :: h, sigma, A
        real(dp), intent(out), dimension(numnod)        :: p0
        END SUBROUTINE
END INTERFACE

!*************************************************************************
!*************************************************************************

REAL(DP)             :: epsilon0                ! to determine support
integer(I4B)         :: i, j, ii, jj, kk        ! loop indices

! Problem domain in 2-dimension:
real(DP)             :: xlowlim, xuplim, ylowlim, yuplim    ! x and y limits
real(DP)             :: L, H                                ! dimesnison of doamin
real(dp)             :: xspac, yspac            ! x & y direction nodal spacing
real(dp)              :: Ra, beta               ! Support radius
real(dp), allocatable :: nodcoord(:,:), nodes(:,:)   ! nodal coordinates (matrix)
real(dp), allocatable :: interval(:,:)
integer(I4B), allocatable  :: conn(:,:)         ! nodal connectivity (matrix)
                                                
! defining a finer grid for plotting:                                                
integer(I4B)  :: nfinedivl, nfinedivw           ! # of mesh points for plotting
integer(I4B)  :: pointnum                       ! point numbers in this mesh
real(dp)      :: smallh, hfinex, hfiney         ! Charactersitic nodal spacing

! Numerical integration parameters:
integer(I4B)  :: quado                         ! order of integartion
integer(I4B)  :: numnod, ntriangle, numq2, numnghbr, nquad    
                                                ! # of nodes, cells, # one/two                    
                                                ! dimensional quadrature points
real(dp), allocatable :: gauss(:,:), pgauss(:,:), wgauss(:)                                                         
real(dp), allocatable :: gs(:,:)
real(dp), dimension(2):: gpos
real(dp)              :: jac, weight

! FPK Equation specific properties:
! Oscillator properties:
real(dp), allocatable         :: diffusion(:,:)   ! This can be a number as well (isotropic case)
real(dp), allocatable         :: drift(:)
real(dp)                      :: zeta, omegan, diff_val

! parameters used by LAPACK equation solver:                                                
integer(I4B)  :: ldKglobal, ldf, info           ! Leading dimension of Kglobal and f
integer(I4B), allocatable, dimension(:) :: v, pivot
real(dp), allocatable              :: Xa(:,:)
integer(I4B)                       :: maxit, ok, refinement
real(dp)                           :: eps
integer(I4B)                       :: nsdim=2                  ! number of space-dimension
CHARACTER(80)                      :: weightfuncname = "gaussian-rbf"
CHARACTER(80)                      :: priorweight
CHARACTER(80)                      :: objectivefuncname = "jaynes"
real(dp), allocatable              :: phi(:)                 ! basis functions
real(dp), allocatable              :: dphi(:,:)              ! first derivatives of basis functions
real(dp), allocatable              :: klocal(:,:), klocalI(:,:), klocalII(:,:)
real(dp), allocatable              :: klocalIII(:,:), klocalIV(:,:)
integer(I4B)                       :: temprow, tempcol
real(dp), allocatable              :: Kglobal(:,:)
real(dp), allocatable              :: Mglobal(:,:), mlocal(:,:)
integer(I4B)                       :: inddirich, lthu
integer(I4B), allocatable          :: essbcnod(:)
real(dp), allocatable              :: ubar(:), f(:), fglobal(:)
real(dp)                           :: q
real(dp), allocatable              :: xfine(:), yfine(:)
real(dp)                           :: xlwlm
real(dp)                           :: xyfine(2)
real(dp), allocatable              :: nodalparam(:)
real(dp)                           :: patxy, p0atxy, Uh, ddx_Uh, ddy_Uh
real(dp), allocatable              :: Uh_gp(:)
integer(I4B)                       :: sizeselectnod
integer(I4B), allocatable          :: selectnod(:)
integer(i4b)                       :: numpoints
real(dp), allocatable              :: Uanlytcl(:)           ! Analytical Solution
real(dp), allocatable              :: grad_Uanlytcl(:,:)    ! Gradient
real(dp)                           :: eh1 = 0.d0            ! H1-semi-norm
real(dp)                           :: el2 = 0.d0            ! L2 norm
real(dp)                           :: eh1_rel = 0.d0        ! Relative H1-semi-norm
real(dp)                           :: el2_rel = 0.d0        ! Relative L2 norm
real(dp)                           :: uh1 = 0.d0        ! H1-semi-norm of analytical
real(dp)                           :: ul2 = 0.d0        ! L2 norm of analytical

! Transient analysis set-up:
real(dp)                           :: dt
real(dp)                           :: A = 5.d0
real(dp)                           :: sigma = 1.d0/9.d0
real(dp)                           :: alpha = 0.5d0         ! Trapezoidal/Cranck Nicolson
integer(i4b)                       :: NTIMESTEPS
real(dp), allocatable              :: p0(:)
real(dp), allocatable              :: v0(:)
real(dp), allocatable              :: p1(:)
real(dp), allocatable              :: v1(:)
real(dp), allocatable              :: p_predictor(:)
real(dp), allocatable              :: Kpred(:)
real(dp), allocatable              :: Kp0(:)
real(dp), allocatable              :: invMglobal(:,:)
real(dp), allocatable              :: M_mod(:,:)
real(dp), allocatable              :: invM_mod(:,:)
integer(i4b)                       :: deltanode     ! node # where delta loading is applied
real(dp), allocatable              :: x0(:)         ! coordinate of node = deltanode

! Various Files:
! Storing the input load, BC, loading, and solution:
character(80)                      :: outfilename = "FPKE2Dmaxent.dat"
! Storing input data (nodal coord, Gauss point, support) to plot in matlab
character(80)                      :: matfilename
! Storing the solution on the finer mesh points at t = 0:
character(80)                      :: solnplott0 = "solnt0.dat"
! Storing the solution on the finer mesh points at t = t_n:
character(80)                      :: solnplottn = "solntn.dat"
! Storing the basis function values at selected nodes:
character(80)                      :: plotbasis = "phi.dat"
! Storing the basis function values at selected nodes:
character(80)                      :: plotgradbasis = "dphi.dat"
! Mesh and EBC data files:
character(80)                      :: meshfile, ebcfile

!*************************************************************************
!*************************************************************************

! Main program begins:
! Defining the Problem Domain:
! Please refer to: http://www.inf.ethz.ch/personal/tulink/FEM14/Ch1_ElmanSyvesterWathen_Ox05.pdf
! Considering simple bi-unit square domain $\Omega$. This is used so as to compare with
! the available analytical solution (inifnite series solution) in the above reference.
! Better to have an external meshfile ready for any dimension nsdim and any number of nodes.

! Domain:
! Read Mesh data from file:
write(*,*)""
write(*,*)"Begining Execution..."
write(*,*)""
write(*,'(A)',ADVANCE="NO") "Please specify the Mesh-data file name: "
read(*,*) meshfile
OPEN(UNIT=51,file = meshfile,STATUS='OLD',ACTION='READ')   ! read from meshfile
read(51,*,END=200)                    ! header: NSD and NODES
read(51,*,END=200)nsdim,numnod        ! dimension and # of nodes
read(51,*,END=200)                    ! header: # of elements
read(51,*,END=200)ntriangle           ! number of "triangular elements"
allocate ( nodcoord(nsdim,numnod) )
allocate ( nodes(nsdim+1,numnod) )
allocate ( conn(3,ntriangle) )
read(51,*,END=200)                    ! header: COORDS
! Nodal Coordinates:
do i = 1,numnod
  read(51,*,END=200)nodes(:,i)
enddo
nodcoord = nodes(2:3,:)
! Nodal Connectivity:
read(51,*,END=200)                    ! header: CONNECTIVITY
do i = 1,ntriangle
  read(51,*,END=200)conn(:,i)
enddo
read(51,*,END=200)              ! header: MAXITER AND TOL
read(51,*,END=200)maxit,eps
read(51,*,END=200)              ! header: PRIOR AND OBJECTIVE
read(51,*)weightfuncname, objectivefuncname
200 CLOSE(UNIT=51)
write(*,*)""
write(*,*)"Read mesh file successfully!"
write(*,*)" "

! Domain:
xlowlim = minval(nodcoord(1,:))
xuplim = maxval(nodcoord(1,:))
ylowlim = minval(nodcoord(2,:))
yuplim = maxval(nodcoord(2,:))
L = xuplim - xlowlim
H = yuplim - ylowlim

! Spatial resolution:
allocate ( interval(numnod,numnod) )
do i=1,numnod
    do j=1,size(nodcoord,2)         ! == numnod
        interval(i,j) = sqrt((nodcoord(1,i)-nodcoord(1,j))**2 + &
                            (nodcoord(2,i)-nodcoord(2,j))**2)
        if (interval(i,j) == 0.d0) then
            interval(i,j) = sqrt(L**2 + H**2)
        endif
    enddo
enddo

h = sqrt(2.d0)*minval(interval)
Ra = 3.5d0*h

! Truncation parameter for Gaussian-rbf prior:
epsilon0 = 1.d-08

! Setting up Quadrature Triangles:
write(*,'(A)',ADVANCE="NO")"Specify Order of Accuracy (can be- 0, 2, 4, 7, 8, or 10):  "
read(*,*) nquad
write(*,*)""

select case (nquad)
    case(0)
        quado = 1
    case(2)
        quado = 3
    case(4)
        quado = 6
    case(7)
        quado = 13
    case(8)
        quado = 16
    case(10)
        quado = 25
    case default
        write(*,*)"ERROR!! NOT CODED YET!"
        STOP
end select

allocate ( gauss (3, quado) )
CALL gaussq_triangle(quado,gauss)

! NOTE: This is not tensor-product quadrature rule
numq2 = ntriangle*quado         ! Total # of quadrature points in 2-dimensions 
allocate ( gs(4,numq2) )       

! Sets up the coordinates for Gauss points in physical space:
CALL egauss(gs,nsdim,nodcoord,conn,gauss,ntriangle,numnod,numq2,quado)

! Oscillator properties:
write(*,'(A)',ADVANCE="NO") "Please specify the Diffusion Coefficient Value: "
read(*,*) diff_val
write(*,*)" "

allocate ( diffusion(nsdim,nsdim) )
!diffusion = 2.d0*diff_val
diffusion = 0.d0
diffusion(nsdim,nsdim) = 2.d0*diff_val

write(*,'(A)',ADVANCE="NO") "Please specify the damping ratio: "
read(*,*)zeta
write(*,*)" "

write(*,*)" "
write(*,'(A)',ADVANCE="NO") "Please specify the natural frequency: "
read(*,*)omegan
write(*,*)" "

priorweight = weightfuncname
write(*,*)""
write(*,'(A)',ADVANCE="NO") "Please specify the EBC data file: "
read(*,*)ebcfile
OPEN(UNIT=52,file=ebcfile,STATUS='OLD',ACTION='READ')   ! read ebc file
read(52,*,END=210)                    ! header: # of EBC nodes
read(52,*,END=210)lthu
allocate ( essbcnod(lthu) )
allocate ( ubar(lthu) )
read(52,*,END=210)                    ! header: NODES & VAL
do i = 1,lthu
  read(52,*,END=210)essbcnod(i), ubar(i) 
enddo
210 CLOSE(UNIT=52)
write(*,*)" "
write(*,*)"Read EBC data file successfully!"
write(*,*)""


write(*,'(A)',ADVANCE="NO")"Specify the Node number to apply the Dirac-delta IC: "
read(*,*) deltanode
write(*,*)" "

allocate ( x0(nsdim) )
do i =1,numnod
    if (nodes(1,i) == deltanode) then
        x0 = nodcoord(:,i)
    endif
enddo
write(*,*)""
write(*,'(A)',ADVANCE="NO")"Specify the time-step (dt) (0.001d0 is good enough): "
read(*,*) dt
write(*,*)""

write(*,*)""
write(*,'(A)',ADVANCE="NO")"Specify the # of time-steps to solve: "
read(*,*) NTIMESTEPS
write(*,*)""

write(*,*)" "
write(*,*) "***************************************************************************************"
write(*,*)" "
write(*,*)" "
write(*,*) "The BVP is set up to solve transient Fokker-Planck-Kolmogorov equation. The domain"
write(*,*) "under consideration is [-10, 10] X [-10, 10] (unstructured grid). The equation is set-up" 
write(*,*) "for 1-dimensional linear oscillator with Gaussian white noise as input."
write(*,*)" "
write(*,*)" "
write(*,*) "***************************************************************************************"

! Computing the Global Stiffness Matrix:
! Initializing Global stiffness matrix: Kglobal
allocate ( Kglobal(numnod,numnod) )
Kglobal = 0.d0

allocate ( Mglobal(numnod,numnod) )
Mglobal = 0.d0

! Setting up some load to solve the system of equation:
allocate ( f(numnod) )
f = 0.d0
q = 0.d0
write(*,*)" "
write(*,*) "Computing the Global Stiffness & Mass Matrix..."
write(*,*) "..."
write(*,*) " "

! Loop over the Gauss points:
do i = 1,size(gs,2)                           ! runs through each column of gs
    gpos = gs(1:2,i)
    weight = gs(3,i)
    jac = gs(4,i)

    ! The drift vector
    allocate ( drift(nsdim) )
    drift(1) = gpos(2)
    drift(2) = -omegan**2*gpos(1) - 2.d0*zeta*omegan*gpos(2)

    ! Determine Nodes in Neighborhood of each Gauss Point (gpos)
    call domain(v,gpos,nodcoord,Ra,numnod)
    
    ! # of neighbors:
    numnghbr = size(v)

    if (numnghbr .lt. 3) then
        write(2,*)"Too few points in the neighborhood of Gauss point # ", i, " a t", gpos, "."
        write(2,*)"Please change the parameters to increase the support size."
        write(6,*)"Too few points in the neighborhood of Gauss point # ", i, " a t", gpos, "."
        write(6,*)"Please change the parameters to increase the support size."
        STOP
    endif

    write(2,*) "Gauss Point #", i, "with coordinates", gpos, "has", numnghbr, "neighbors."
    write(2,*) "Their node #'s are: (", (v(kk), kk=1,numnghbr), ")"

    allocate ( Xa(numnghbr,nsdim) )
    call getcoordinates(Xa,v,numnghbr,numnod,nodcoord)

    allocate ( phi(numnghbr) )
    allocate ( dphi(numnghbr,nsdim) )
    
    CALL driver(numnghbr, nsdim, weightfuncname, epsilon0, objectivefuncname, &
                xa, gpos, Ra, maxit, eps, phi, dphi)
        
    ! Initializing "Local" stiffness matrix: klocal. Size of this varies for each Gauss point
    allocate ( klocal(numnghbr,numnghbr) )
    allocate ( mlocal(numnghbr,numnghbr) )

    allocate ( klocalI(nsdim,numnghbr) )

    klocalI = matmul(diffusion, transpose(dphi))

    allocate ( klocalII(numnghbr,numnghbr) )
    klocalII = 0.5d0*matmul(dphi, klocalI)      ! Goes as the second term of stiffness matrix

    allocate ( klocalIII(nsdim, numnghbr) )
    
    do ii=1,nsdim
        do jj=1,numnghbr
            klocalIII(ii,jj) = drift(ii)*phi(jj)     ! b*transpose(phi)
        enddo
    enddo

    allocate ( klocalIV(numnghbr,numnghbr) )
    klocalIV = matmul(dphi, klocalIII)              ! Goes as the second term of stiffness matrix

    
    klocal = 0.5d0*weight*jac*(klocalII - klocalIV)    ! local stiffness: value of the integrand
                                                            ! evaluated at this Gauss point.
    do ii=1,numnghbr
        do jj =1,numnghbr
            mlocal(ii,jj) = phi(ii)*phi(jj)
        enddo
    enddo
    mlocal = 0.5d0*weight*jac*mlocal

    deallocate (klocalI)
    deallocate (klocalII)
    deallocate (klocalIII)
    deallocate (klocalIV)
    deallocate (drift)

    ! Local to global number mapping:
    ! In local neighbourhood of gpos- node # (again local) ii is mapped to global node # v(ii)
    do ii = 1, numnghbr
        temprow = v(ii)                 ! Global row number corresponding to node # v(ii)
        do jj = 1, numnghbr
            ! Node # v(ii) --> local sequential index ii
            tempcol = v(jj)             ! global column number corresponding to node # v(jj)
            Kglobal(temprow,tempcol) = Kglobal(temprow,tempcol) + klocal(ii,jj)     ! This would be singular
            Mglobal(temprow,tempcol) = Mglobal(temprow,tempcol) + mlocal(ii,jj)
        enddo
        f(temprow) = f(temprow) + 0.5d0*weight*jac*phi(ii)*q
    enddo

    deallocate (phi)
    deallocate (dphi)
    deallocate (Xa)
    deallocate (v)
    deallocate (klocal)
    deallocate (mlocal)
enddo
write(*,*)" "
write(*,*) "Finished Computing the Global Stiffness & Mass Matrix!"
write(*,*) "******************************************************"
write(*,*)" "

! Modifying global stiffness matrix and consistent load vector to incorporate EBC:
do kk = 1,size(essbcnod)
    do ii = 1,size(Kglobal,1)
        if (ii .ne. essbcnod(kk)) then
             f(ii) = f(ii) - Kglobal(ii,essbcnod(kk))*ubar(kk)
         endif
    enddo
    Kglobal(essbcnod(kk),:) = 0.d0
    Mglobal(essbcnod(kk),:) = 0.d0
    Kglobal(:,essbcnod(kk)) = 0.d0
    Mglobal(:,essbcnod(kk)) = 0.d0
    Kglobal(essbcnod(kk),essbcnod(kk)) = 1.d0
    Mglobal(essbcnod(kk),essbcnod(kk)) = 1.d0
    f(essbcnod(kk)) = ubar(kk)
enddo

! Storing the global load vector for output:
allocate ( fglobal(numnod) )
fglobal = f

allocate ( p0(numnod) )
allocate ( v0(numnod) )
allocate ( p1(numnod) )
allocate ( v1(numnod) )
allocate ( p_predictor(numnod) )
allocate ( Kpred(numnod) )
allocate ( Kp0(numnod) )
allocate ( invMglobal(numnod,numnod) )
allocate ( M_mod(numnod,numnod) )
allocate ( invM_mod(numnod,numnod) )

! Defining a very fine mesh:
write(*,*)""
nfinedivl = 1000

hfinex = L/nfinedivl
hfiney = hfinex

xlwlm = xlowlim

! nodal coordinates of this fine mesh (x)
allocate ( xfine(nfinedivl+1) )
do i=1,(nfinedivl+1)
    xfine(i) = xlwlm + (i-1)*hfinex
enddo
yfine = xfine

! Open solution output files:

! Storing the solution on the finer mesh points: "solnt0.dat"
open(UNIT=22, file = solnplott0, STATUS='NEW',iostat=ok)

! Cranck-Nicolson Routine:
CALL initialcond(p0, nodcoord, h, A, sigma, x0, nsdim, numnod)

write(*,*)"Generating MATLAB plot file for solution at time t = 0..."

do j = 1,(nfinedivl + 1)
    do i = 1,(nfinedivl + 1)

        xyfine = (/ xfine(i), yfine(j) /)
        
        ! Determine Nodes in Neighborhood of each Gauss Point (gpos)
        call domain(v,xyfine,nodcoord,Ra,numnod)

        numnghbr = size(v)
        if (numnghbr .lt. 3) then
            write(*,*)"Too few points in the neighborhood of the mesh point"
            stop
        endif
        
        allocate ( Xa(numnghbr,nsdim) )
        call getcoordinates(Xa,v,numnghbr,numnod,nodcoord)

        allocate ( nodalparam(numnghbr) )
        do ii=1,numnghbr
            nodalparam(ii) = p0(v(ii))
        enddo

        allocate ( phi(numnghbr) )
        allocate ( dphi(numnghbr,nsdim) )

        CALL driver(numnghbr, nsdim, weightfuncname, epsilon0, objectivefuncname, &
                    Xa, xyfine, Ra, maxit, eps, phi, dphi)

        p0atxy = dot_product(nodalparam,phi)

        write(22,*) xyfine(1), xyfine(2), p0atxy

        deallocate (phi)
        deallocate (dphi)
        deallocate (Xa)
        deallocate (v)
        deallocate (nodalparam)
    enddo
enddo
write(*,*)"Finished writing the MATLAB plot file for solution at time t = 0."
write(*,*)""

CALL invmatrix(Mglobal,invMglobal,numnod)
Kp0 = matmul(Kglobal, p0)
v0 = matmul(invMglobal, (fglobal - Kp0))

! The following matrices don't change with time:
M_mod = Mglobal + alpha*dt*Kglobal
CALL invmatrix(M_mod, invM_mod, numnod)

write(*,*) "Now solving..."
write(*,*) "..."
write(*,*) " "
do jj=1,NTIMESTEPS      ! Loop over time
    p_predictor = p0 + (1.d0 - alpha)*dt*v0
    Kpred = matmul(Kglobal, p_predictor)
    v1 = matmul(invM_mod, (fglobal - Kpred))
    p1 = p_predictor + alpha*dt*v1
    write(*,*)""
    write(*,*)" Finished Computing solution at time-step = ", jj
    write(*,*)""
    v0 = v1
    p0 = p1
enddo

! Main loop to evaluate the function value at last time point:
! Storing the solution on the finer mesh points: "solntn.dat"
open(UNIT=12, file = solnplottn, STATUS='NEW',iostat=ok)
do j = 1,(nfinedivl + 1)
    do i = 1,(nfinedivl + 1)

        xyfine = (/ xfine(i), yfine(j) /)
        
        ! Determine Nodes in Neighborhood of each Gauss Point (gpos)
        call domain(v,xyfine,nodcoord,Ra,numnod)

        numnghbr = size(v)
        if (numnghbr .lt. 3) then
            write(*,*)"Too few points in the neighborhood of the mesh point"
            stop
        endif
        
        allocate ( Xa(numnghbr,nsdim) )
        call getcoordinates(Xa,v,numnghbr,numnod,nodcoord)

        allocate ( nodalparam(numnghbr) )
        do ii=1,numnghbr
            nodalparam(ii) = p0(v(ii))
        enddo

        allocate ( phi(numnghbr) )
        allocate ( dphi(numnghbr,nsdim) )

        CALL driver(numnghbr, nsdim, weightfuncname, epsilon0, objectivefuncname, &
                    Xa, xyfine, Ra, maxit, eps, phi, dphi)

        patxy = dot_product(nodalparam,phi)

        write(12,*) xyfine(1), xyfine(2), patxy

        deallocate (phi)
        deallocate (dphi)
        deallocate (Xa)
        deallocate (v)
        deallocate (nodalparam)
    enddo
enddo
write(*,*)""
write(*,*) "***************************************************************************************"
write(*,*) "************************************** SOLVED! ****************************************"
write(*,*) "***************************************************************************************"
write(*,*)""

!*************************************************************************
!*************************************************************************

! Writing the ouput to a file:
! opening the output file:
! filename = FPKE2Dmaxent.dat
write(*,*) "Now writing output to an external data file..."
write(*,*) "..."
write(*,*) " "

open(UNIT=2, file =  outfilename, STATUS='NEW')
write(2,*)"Spatial dimension = ", nsdim
write(2,*)"Total number of nodes = ", numnod
write(2,*)"Number of Gauss points in 2D = ", numq2
write(2,*)"Number of integration cells = ", ntriangle
write(2,*)"Nodal Spacing = ", h
write(2,*)" "

! Nodal Coordinates:
write(2,*)"Nodal Coordinates: "
do i = 1,numnod
  write(2,*)"Node = ",i," Coordinates = ", nodcoord(:,i)
enddo
write(2,*)" "

! Dirichlet BC nodes:
write(2,*)"The Nodes with imposed EBC: "
do i = 1, size(essbcnod)
    kk = essbcnod(i)
    write(2,*)"Node = ",kk," Coordinates = ", nodcoord(:,kk)
enddo
write(2,*)" "

! Problem Setup details:
write(2,*)"Maximum number of iterations = ", maxit
write(2,*)"Convergence tolerance = ", eps
write(2,*)"Support radius = ", Ra
write(2,*)"Prior weight function: ", weightfuncname
write(2,*)" "

!*************************************************************************
!*************************************************************************

! Plotting (writing another MATLAB inputfile to plot the solution and 
! few shape functions:

! Storing the basis function values at selected nodes:  "phi.dat"
open(UNIT=16, file = plotbasis, STATUS='NEW',iostat=ok)

! Storing the derivative of basis function values at selected nodes: "dphi.dat"
open(UNIT=17, file = plotgradbasis, STATUS='NEW',iostat=ok)

write(*,*)" "
write(*,'(A)',ADVANCE="NO")"Specify the number of nodes to plot:  "
read(*,*) sizeselectnod
allocate ( selectnod(sizeselectnod) )

write(*,*)" "
write(*,*)"NOTE: Total # of nodes used for computation =  ", numnod
write(*,'(A)',ADVANCE="NO")"Specify the nodes (numbers) to plot: "
read(*,*) selectnod

write(*,*)" "
write(*,*)"Computing the solution over the finer mesh..."
write(*,*)" "

! Main loop to evaluate the function value:
do j = 1,(nfinedivl + 1)
    do i = 1,(nfinedivl + 1)

        xyfine = (/ xfine(i), yfine(j) /)
        
        ! Determine Nodes in Neighborhood of each Gauss Point (gpos)
        call domain(v,xyfine,nodcoord,Ra,numnod)

        numnghbr = size(v)
        if (numnghbr .lt. 3) then
             write(*,*)"Too few points in the neighborhood of the mesh point"
             stop
        endif
        
        allocate ( Xa(numnghbr,nsdim) )
        call getcoordinates(Xa,v,numnghbr,numnod,nodcoord)

        allocate ( nodalparam(numnghbr) )
        do ii=1,numnghbr
            nodalparam(ii) = f(v(ii))
        enddo

        allocate ( phi(numnghbr) )
        allocate ( dphi(numnghbr,nsdim) )

        CALL driver(numnghbr, nsdim, weightfuncname, epsilon0, objectivefuncname, &
                    Xa, xyfine, Ra, maxit, eps, phi, dphi)


        do ii=1,size(v)
            do jj=1,sizeselectnod
                if (v(ii) == selectnod(jj)) then
                    write(16,*) xyfine(1), xyfine(2), phi(ii)
                    write(17,*) xyfine(1), xyfine(2), dphi(ii,:)
                endif
            enddo
        enddo

        deallocate (phi)
        deallocate (dphi)
        deallocate (Xa)
        deallocate (v)
        deallocate (nodalparam)
    enddo
enddo

!*************************************************************************
!*************************************************************************

!write(*,*)"Computing the H1 semi-norm and L2 norm..."
!write(*,*)"..."
!write(*,*)
!The Analytical Solution and its derivatives:
!allocate ( Uanlytcl(numq2) )
!allocate ( grad_Uanlytcl(numq2, nsdim) )
!CALL anlytclsol(Uanlytcl, grad_Uanlytcl, gs, numq2, nsdim)

! Loop over the finer Gauss points:
!do i = 1,size(gs,2)           ! runs through each column of gs
!    gpos = gs(1:2,i)
!    weight = gs(3,i)
!    jac = gs(4,i)
!
!    ! Determine Nodes in Neighborhood of each Gauss Point (gpos)
!    call domain(v, gpos, nodcoord, Ra, numnod)
!    
!    ! # of neighbors:
!    numnghbr = size(v)
!
!    if (numnghbr .lt. 3) then
!        write(*,*)"Too few points in the neighborhood of finer Gauss point # ", i, " at", gpos, "."
!        write(*,*)"Please change the parameters to increase the support size."
!        STOP
!    endif
!
!    allocate ( Xa(numnghbr,nsdim) )
!    call getcoordinates(Xa,v,numnghbr,numnod,nodcoord)
!
!    allocate ( nodalparam(numnghbr) )
!
!    do ii=1,numnghbr
!        nodalparam(ii) = f(v(ii))
!    enddo
!
!    allocate ( phi(numnghbr) )
!    allocate ( dphi(numnghbr,nsdim) )
!
!    CALL driver(numnghbr, nsdim, weightfuncname, epsilon0, objectivefuncname, &
!                    Xa, gpos, Ra, maxit, eps, phi, dphi)
!
!    ! Uh, ddx_Uh, ddy_Uh :
!    Uh = dot_product(nodalparam,phi)
!    ddx_Uh = dot_product(nodalparam,dphi(:,1))
!    ddy_Uh = dot_product(nodalparam,dphi(:,2))
!
!    ! Computing error-norms:
!    ! H1 semi-norm (running sum):
!
!    uh1 = uh1 +  0.5d0*weight*jac*(grad_Uanlytcl(i,1)**2 + grad_Uanlytcl(i,2)**2)
!
!    eh1 = eh1 +  0.5d0*weight*jac*((grad_Uanlytcl(i,1) - ddx_Uh)**2 + &
!                                (grad_Uanlytcl(i,2) - ddy_Uh)**2)
!
!    ! L2-norm (running sum):
!    ul2 = ul2 + weight*jac*Uanlytcl(i)**2
!    el2 = el2 + weight*jac*(Uanlytcl(i) - Uh)**2
!    
!    deallocate (nodalparam)
!    deallocate (phi)
!    deallocate (dphi)
!    deallocate (Xa)
!    deallocate (v)
!enddo

!eh1 = sqrt(eh1)
!uh1 = sqrt(uh1)
!el2 = sqrt(el2)
!ul2 = sqrt(ul2)
!eh1_rel = 100.d0*eh1/uh1
!el2_rel = 100.d0*el2/ul2


!write(*,*) " "
!write(*,*) "a-priori error estimates..."
!write(*,*) "H1-seminorm of error in solution = ", eh1
!write(*,*) "H1-seminorm of relative error in solution = ", eh1_rel, "%"

!write(*,*) " "
!write(*,*) "L2norm of error in solution = ", el2
!write(*,*) "L2norm of relative error in solution = ", el2_rel, "%"
!write(*,*) " "

!write(2,*)"H1-seminorm (AKA: energy semi-norm) = ", eh1
!write(2,*)"H1-seminorm of relative error in solution = ", eh1_rel, "%"
!write(2,*)" "

!write(2,*)"L2norm of error = ", el2
!write(2,*)"L2norm of relative error in solution = ", el2_rel, "%"
!write(2,*)" "


close(UNIT=2)
close(UNIT=22)
close(UNIT=12)
close(UNIT=16)
close(UNIT=17)

write(*,*) "***************************************************************************************"
write(*,*) "**************************************** BYE! *****************************************"
write(*,*) "***************************************************************************************"

END PROGRAM maxent2dFPK

!*************************************************************************
!*************************************************************************

!*************************************************************************
!*************************************************************************

SUBROUTINE  egauss(gs,nsdim,nodcoord,conn,gauss,ntriangle,numnod,numq2,quado)
USE NRTYPE
! function to set up Gauss points, Jacobian, and weights
IMPLICIT NONE
! Interfacing part of the dummy variables
integer(i4b), intent(in)                         :: ntriangle, numnod
integer(i4b), intent(in)                         :: numq2, quado, nsdim
integer(i4b), dimension(3,ntriangle), intent(in) :: conn          ! DELAUNAY CONNECTIVITY (CCW)
real(dp), dimension(3,quado), intent(in)         :: gauss         ! 
real(dp), dimension(nsdim,numnod), intent(in)    :: nodcoord      ! Nodal Coordinate array
real(dp), dimension(4,numq2), intent(out)        :: gs            ! Global coordinates of Gauss points

! Local variables:
integer                     :: index = 0
integer(i4b)                :: l, e, i, j, je
real(dp), allocatable       :: xe(:), ye(:)
real(dp)                    :: eta, xi, jcob, xq, yq, N1, N2, N3

! Nodal coordinates (x,y) of each triangle:
allocate ( xe(3) )
allocate ( ye(3) )

l = size(gauss,2)       ! # of Gauss points for each triangle

do e = 1,ntriangle
    ! Determine Nodes & Nodal Coordinates in Each Triangle
    do  j = 1,3
        je = conn(j,e)      ! node numbers for element # e
        ! nodal coordinate of node "je" or element "e":
        xe(j) = nodcoord(1,je)
        ye(j) = nodcoord(2,je)
    enddo
    
    do i = 1,l
        index = index + 1
        xi = gauss(1,i)
        eta = gauss(2,i)
        ! Triangle Isoparametric Shape functions:
        N1 = 1.d0 - xi - eta
        N2 = xi
        N3 = eta
                       
        xq = N1*xe(1) + N2*xe(2) + N3*xe(3)
        yq = N1*ye(1) + N2*ye(2) + N3*ye(3)

        jcob = abs(xe(1)*(ye(2)-ye(3)) + xe(2)*(ye(3)-ye(1)) &
                        + xe(3)*(ye(1)-ye(2)))      ! 2*area

        gs(1,index) = xq
        gs(2,index) = yq
        gs(3,index) = gauss(3,i)
        gs(4,index) = jcob
    enddo
enddo
END SUBROUTINE
!*************************************************************************
!*************************************************************************

SUBROUTINE domain(v,gpos,x,Ra,nnodes)
USE NRTYPE
! Determines Nodes whose Domain of Influence includes the Point- Gpos (Gauss Point)
implicit none
! Interfacing part of the dummy variables
real(DP), dimension(2), intent(in)                      :: gpos
integer(I4B), intent(in)                                :: nnodes
real(DP), intent(in)                                    :: Ra
real(DP), dimension(2,nnodes), intent(in)               :: x
integer(I4B), allocatable, intent(out), dimension(:)    :: v

! Local variables:
real(DP)                           :: macheps
integer(I4B)                       :: i, counter
real(dp), dimension(nnodes)        :: norm
macheps = epsilon(0d0)                            ! machine precision

do i=1,nnodes
    norm(i) = sqrt( (x(1,i)-gpos(1))**2 + (x(2,i)-gpos(2))**2 )
enddo

counter = 0
do i=1,nnodes
    if (norm(i) .le. Ra) then
        counter = counter + 1
    endif
enddo

allocate ( v(counter) )
counter = 0
do i = 1, nnodes             ! across all nodes/columns/1 to nnodes
    if (norm(i) .le. Ra) then
        counter = counter + 1
        v(counter) = i
    endif
enddo

END SUBROUTINE

!*************************************************************************
!*************************************************************************

SUBROUTINE getcoordinates(xyz,v,numnghbr,nnodes,x)
USE NRTYPE
implicit none

! Interfacing part of the dummy variables
integer(I4B), dimension(numnghbr),  intent(in)      :: v
integer(I4B), intent(in)                            :: nnodes, numnghbr
real(DP), dimension(2,nnodes), intent(in)           :: x
real(DP), dimension(numnghbr,2), intent(out)        :: xyz

! Local variables:
integer(I4B)                                        :: i, j

do i=1,numnghbr
    do j=1,nnodes
        if (v(i) == j) then
            xyz(i,1) = x(1,j)
            xyz(i,2) = x(2,j)
        endif
    enddo
enddo


END SUBROUTINE

!*************************************************************************
!*************************************************************************

SUBROUTINE supportplot(gs,nodcoord,Ra,numq2,nnodes,filename)
USE NRTYPE
implicit none
! Interfacing part of the dummy variables
integer(I4B), intent(in)                     :: numq2, nnodes
real(DP), dimension(4,numq2), intent(out)    :: gs
real(DP), dimension(2,nnodes), intent(in)    :: nodcoord
real(DP), intent(in)                         :: Ra
character(80),  intent(in)                   :: filename

! Local variables:
integer(I4B)                                 :: i, ok

! Writes output file to be read by MATLAB:
! open output file 
open(UNIT=3, file = filename, STATUS='NEW',iostat=ok)
! Total number of nodes:
write(3,*) nnodes
! Total number of quadrature points in 2-dimensions
write(3,*) numq2
! Support radius at each node:
write(3,*) Ra
! Nodal coordinates
do i = 1,nnodes
  write(3,*) i, nodcoord(:,i)
enddo
! Quadrature point coordintaes
do i = 1,numq2
    write(3,*) i, gs(1,i), gs(2,i)
enddo

close(UNIT=3)

END SUBROUTINE

!*************************************************************************
!*************************************************************************

SUBROUTINE gridcircle(allnodes,nodloc,Ra,ndivtheta,ndivr)
USE NRTYPE
implicit none
! Interfacing part of the dummy variables
integer(I4B), intent(in)                                      :: ndivtheta, ndivr
real(DP), intent(in)                                          :: Ra
reaL(DP), dimension(2), intent(in)                            :: nodloc
real(DP), dimension(2,(1 + (ndivtheta+1)*ndivr)), intent(out) :: allnodes

! Local variables:
integer(I4B)                                 :: i, j, nodnumber
real(DP)                                     :: rinc, thetainc

rinc = Ra/ndivr
thetainc = Pi/ndivtheta
do i=1,(ndivtheta+1)
    do j=1,ndivr
        nodnumber = (i-1)*ndivr + j
        allnodes(1,nodnumber) = nodloc(1) + j*rinc*cos((i-1)*thetainc)
        allnodes(2,nodnumber) = nodloc(2) + j*rinc*sin((i-1)*thetainc)
    enddo
enddo
allnodes(:,(1 + (ndivtheta+1)*ndivr)) = nodloc

END SUBROUTINE

!*************************************************************************
!*************************************************************************

SUBROUTINE anlytclsol(Uanlytcl, grad_Uanlytcl, gparray, numpoints, nsdim)

USE NRTYPE
IMPLICIT NONE
! Interfacing part of the dummy variables
integer(i4b), intent(in)                            :: numpoints, nsdim
real(dp), intent(in), dimension(4,numpoints)        :: gparray
real(dp), dimension(numpoints), intent(out)         :: Uanlytcl
real(dp), dimension(numpoints,nsdim), intent(out)   :: grad_Uanlytcl

! Local variables:
integer(i4b)       :: i
real(dp)           :: x, y

do i=1,numpoints
    x = gparray(1,i)
    y = gparray(2,i)
    Uanlytcl(i) = sin(PI*x)*sin(PI*y)
    grad_Uanlytcl(i,1) = PI*cos(PI*x)*sin(PI*y)
    grad_Uanlytcl(i,2) = PI*sin(PI*x)*cos(PI*y)
enddo

END SUBROUTINE anlytclsol

!*************************************************************************
!*************************************************************************
!
!     subroutine gaussq_triangle(nint,pgauss,wgauss)
!     Purpose
!     =======
!     Gauss quadrature points and weights for triangular element
!     input: nint        - number of quadrature points
!     output: pgauss(3,*) - natural coordinates of Gauss points
!             wgauss(*)   - weights of Gauss points
!
!     Comments
!     ========
!     Rules for nint = 1,3,6,13 are taken from the finite
!     element book: Cook, Malkus, Plesha, "Concepts and Applns
!     of the Finite Element Method"
!     Rules for nint >= 20 are from the paper: Dunavant, D. A.,
!     "High Degree Efficient Symmetrical Gaussian Quadrature
!     Rules for the Triangle", IJNME, vol. 21, 1129-1148, 1985
!
!***************************************************************
!
SUBROUTINE gaussq_triangle(nint,gauss)
USE NRTYPE
IMPLICIT NONE
! Interfacing part of the dummy variables
INTEGER(I4B), INTENT(IN)    :: nint
!  1st row gives the xi-coordinates of points %
!  2nd row gives the eta-coordinates of points %
!  3rd row gives the weights
REAL(DP), INTENT(OUT)       :: gauss(3,nint)

! Local Variable:
integer i

if (nint .eq. 1) then
! zeroth-order accurate for triangular element
    gauss(1,1) = 0.33333333333333d0
    gauss(2,1) = 0.33333333333333d0
    gauss(3,1) = 1.00000000000000d0

elseif (nint .eq. 3) then
! second-order accurate for triangular element
    gauss(1,1) = 0.16666666666667d0
    gauss(2,1) = 0.16666666666667d0
    gauss(3,1) = 0.33333333333333d0

    gauss(1,2) = 0.16666666666667d0
    gauss(2,2) = 0.66666666666667d0
    gauss(3,2) = 0.33333333333333d0

    gauss(1,3) = 0.66666666666667d0
    gauss(2,3) = 0.16666666666667d0
    gauss(3,3) = 0.33333333333333d0

elseif (nint .eq. 6) then
! fourth-order accurate for triangular element
    gauss(1,1) = 0.44594849091597d0  
    gauss(2,1) = 0.44594849091597d0
    gauss(3,1) = 0.22338158967801d0

    gauss(1,2) = 0.44594849091597d0  
    gauss(2,2) = 0.10810301816807d0
    gauss(3,2) = 0.22338158967801d0

    gauss(1,3) = 0.10810301816807d0  
    gauss(2,3) = 0.44594849091597d0
    gauss(3,3) = 0.22338158967801d0

    gauss(1,4) = 0.09157621350977d0  
    gauss(2,4) = 0.09157621350977d0
    gauss(3,4) = 0.10995174365532d0

    gauss(1,5) = 0.09157621350977d0 
    gauss(2,5) = 0.81684757298046d0
    gauss(3,5) = 0.10995174365532d0

    gauss(1,6) = 0.81684757298046d0
    gauss(2,6) = 0.09157621350977d0
    gauss(3,6) = 0.10995174365532d0

elseif (nint .eq. 13) then
! seventh-order accurate for triangular element
    gauss(1,1) = 0.33333333333333d0
    gauss(2,1) = 0.33333333333333d0
    gauss(3,1) = -0.14957004446768d0

    gauss(1,2) = 0.26034596607904d0  
    gauss(2,2) = 0.26034596607904d0
    gauss(3,2) = 0.17561525743321d0

    gauss(1,3) = 0.26034596607904d0
    gauss(2,3) = 0.47930806784192d0
    gauss(3,3) = 0.17561525743321d0

    gauss(1,4) = 0.47930806784192d0
    gauss(2,4) = 0.26034596607904d0
    gauss(3,4) = 0.17561525743321d0

    gauss(1,5) = 0.06513010290222d0  
    gauss(2,5) = 0.06513010290222d0
    gauss(3,5) = 0.05334723560884d0

    gauss(1,6) = 0.06513010290222d0  
    gauss(2,6) = 0.86973979419557d0
    gauss(3,6) = 0.05334723560884d0

    gauss(1,7) = 0.86973979419557d0  
    gauss(2,7) = 0.06513010290222d0
    gauss(3,7) = 0.05334723560884d0

    gauss(1,8) = 0.31286549600487d0  
    gauss(2,8) = 0.63844418856981d0
    gauss(3,8) = 0.07711376089026d0

    gauss(1,9) = 0.63844418856981d0  
    gauss(2,9) = 0.04869031542532d0
    gauss(3,9) = 0.07711376089026d0

    gauss(1,10) = 0.04869031542532d0  
    gauss(2,10) = 0.31286549600487d0
    gauss(3,10) = 0.07711376089026d0

    gauss(1,11) = 0.63844418856981d0 
    gauss(2,11) = 0.31286549600487d0
    gauss(3,11) = 0.07711376089026d0

    gauss(1,12) = 0.31286549600487d0
    gauss(2,12) = 0.04869031542532d0
    gauss(3,12) = 0.07711376089026d0

    gauss(1,13) = 0.04869031542532d0
    gauss(2,13) = 0.63844418856981d0
    gauss(3,13) = 0.07711376089026d0

elseif (nint .eq. 16) then
! eighth-order accurate for triangular element
    gauss(1,1) = 0.33333333333333d0
    gauss(2,1) = 0.33333333333333d0
    gauss(3,1) = 0.14431560767779d0

    gauss(1,2) = 0.45929258829272d0
    gauss(2,2) = 0.45929258829272d0
    gauss(3,2) = 0.09509163426728d0

    gauss(1,3) = 0.45929258829272d0 
    gauss(2,3) = 0.08141482341455d0
    gauss(3,3) = 0.09509163426728d0

    gauss(1,4) = 0.08141482341455d0
    gauss(2,4) = 0.45929258829272d0
    gauss(3,4) = 0.09509163426728d0

    gauss(1,5) = 0.17056930775176d0 
    gauss(2,5) = 0.17056930775176d0
    gauss(3,5) = 0.10321737053472d0

    gauss(1,6) = 0.17056930775176d0
    gauss(2,6) = 0.65886138449648d0
    gauss(3,6) = 0.10321737053472d0

    gauss(1,7) = 0.65886138449648d0  
    gauss(2,7) = 0.17056930775176d0
    gauss(3,7) = 0.10321737053472d0

    gauss(1,8) = 0.05054722831703d0
    gauss(2,8) = 0.05054722831703d0
    gauss(3,8) = 0.03245849762320d0

    gauss(1,9) = 0.05054722831703d0
    gauss(2,9) = 0.89890554336594d0
    gauss(3,9) = 0.03245849762320d0

    gauss(1,10) = 0.89890554336594d0
    gauss(2,10) = 0.05054722831703d0
    gauss(3,10) = 0.03245849762320d0

    gauss(1,11) = 0.26311282963464d0
    gauss(2,11) = 0.72849239295540d0
    gauss(3,11) = 0.02723031417443d0

    gauss(1,12) = 0.72849239295540d0
    gauss(2,12) = 0.00839477740996d0
    gauss(3,12) = 0.02723031417443d0

    gauss(1,13) = 0.00839477740996d0
    gauss(2,13) = 0.26311282963464d0
    gauss(3,13) = 0.02723031417443d0

    gauss(1,14) = 0.72849239295540d0
    gauss(2,14) = 0.26311282963464d0
    gauss(3,14) = 0.02723031417443d0

    gauss(1,15) = 0.26311282963464d0  
    gauss(2,15) = 0.00839477740996d0
    gauss(3,15) = 0.02723031417443d0

    gauss(1,16) = 0.00839477740996d0  
    gauss(2,16) = 0.72849239295540d0
    gauss(3,16) = 0.02723031417443d0
elseif (nint .eq. 25) then
! tenth-order accurate for triangular element
    gauss(1,1) = 1.d0/3.d0
    gauss(2,1) = gauss(1,1)
    gauss(3,1) = 0.090817990382754d0

    gauss(1,2) = 0.028844733232685d0
    gauss(2,2) = 0.485577633383657d0
    gauss(3,2) = 0.036725957756467d0

    gauss(1,3) = gauss(2,2)
    gauss(2,3) = gauss(1,2)
    gauss(3,3) = 0.036725957756467d0

    gauss(1,4) = gauss(2,2)
    gauss(2,4) = gauss(2,2)
    gauss(3,4) = 0.036725957756467d0

    gauss(1,5) = 0.781036849029926d0
    gauss(2,5) = 0.109481575485037d0
    gauss(3,5) = 0.045321059435528d0

    gauss(1,6) = gauss(2,5)
    gauss(2,6) = gauss(1,5)
    gauss(3,6) = 0.045321059435528d0

    gauss(1,7) = gauss(2,5)
    gauss(2,7) = gauss(2,5)
    gauss(3,7) = 0.045321059435528d0

    gauss(1,8) = 0.141707219414880d0
    gauss(2,8) = 0.307939838764121d0
    gauss(3,8) = 0.072757916845420d0

    gauss(1,9) = gauss(1,8)
    gauss(2,9) = gauss(3,8)
    gauss(3,9) = 0.072757916845420d0

    gauss(1,10) = gauss(2,8)
    gauss(2,10) = gauss(1,8)
    gauss(3,10) = 0.072757916845420d0

    gauss(1,11) = gauss(2,8)
    gauss(2,11) = gauss(3,8)
    gauss(3,11) = 0.072757916845420d0
    
    gauss(1,12) = gauss(3,8)
    gauss(2,12) = gauss(1,8)
    gauss(3,12) = 0.072757916845420d0

    gauss(1,13) = gauss(3,8)
    gauss(2,13) = gauss(2,8)
    gauss(3,13) = 0.072757916845420d0

    gauss(1,14) = 0.025003534762686d0
    gauss(2,14) = 0.246672560639903d0
    gauss(3,14) = 0.028327242531057d0

    gauss(1,15) = gauss(1,14)
    gauss(2,15) = gauss(3,14)
    gauss(3,15) = 0.028327242531057d0

    gauss(1,16) = gauss(2,14)
    gauss(2,16) = gauss(1,14)
    gauss(3,16) = 0.028327242531057d0

    gauss(1,17) = gauss(2,14)
    gauss(2,17) = gauss(3,14)
    gauss(3,17) = 0.028327242531057d0

    gauss(1,18) = gauss(3,14)
    gauss(2,18) = gauss(1,14)
    gauss(3,18) = 0.028327242531057d0

    gauss(1,19) = gauss(3,14)
    gauss(2,19) = gauss(2,14)
    gauss(3,19) = 0.028327242531057d0

    gauss(1,20) = 0.009540815400299d0
    gauss(2,20) = 0.066803251012200d0
    gauss(3,20) = 0.009421666963733d0

    gauss(1,21) = gauss(1,20)
    gauss(2,21) = gauss(3,20)
    gauss(3,21) = 0.009421666963733d0

    gauss(1,22) = gauss(2,20)
    gauss(2,22) = gauss(1,20)
    gauss(3,22) = 0.009421666963733d0

    gauss(1,23) = gauss(2,20)
    gauss(2,23) = gauss(3,20)
    gauss(3,23) = 0.009421666963733d0

    gauss(1,24) = gauss(3,20)
    gauss(2,24) = gauss(1,20)
    gauss(3,24) = 0.009421666963733d0

    gauss(1,25) = gauss(3,20)
    gauss(2,25) = gauss(2,20)
    gauss(3,25) = 0.009421666963733d0

endif
END SUBROUTINE gaussq_triangle

!*************************************************************************
!*************************************************************************

SUBROUTINE initialcond(p0, nodcoord, h, A, sigma, x0, nsdim, numnod)
USE NRTYPE
IMPLICIT NONE
! Interfacing part of the dummy variables
integer(i4b), intent(in)                        :: nsdim, numnod
real(dp), intent(in), dimension(nsdim,numnod)   :: nodcoord
real(dp), intent(in), dimension(nsdim)          :: x0
real(dp), intent(in)                            :: h, sigma, A
real(dp), intent(out), dimension(numnod)        :: p0

! Local variables:
integer(i4b)       :: i
real(dp)           :: r, fval

do i=1,numnod
    r = (0.5d0*(nodcoord(1,i)-x0(1))**2/sigma**2 + & 
        0.5d0*(nodcoord(2,i)-x0(2))**2/sigma**2)
    fval = A*exp(-r)
    if (fval <= 1.d-06) then
        p0(i) = 0.d0
    else 
        p0(i) = fval
    endif
enddo

END SUBROUTINE initialcond
