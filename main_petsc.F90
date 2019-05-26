program main_petsc

  use fd

  implicit none  

#include <finclude/petscdef.h>
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
  !!  use petsc
  !!  use mpi

  type(mesh2d) :: Th
  type(dd1d)   :: ddgeo
  type(OSMopt) :: OSM

  double precision, allocatable :: fullAsp(:,:)
  double precision, allocatable :: vecX(:), vecY(:)
  
  type(spmat)  :: Asp
  type(spmat)  :: AaugD, AaugOff

  type(spmat), allocatable  :: res(:), resGAM(:)
  type(spmat), allocatable  :: Asub(:), Asub2intf(:,:), Aintf(:)

  type(sp_csr) :: Bsp
  type(sp_csr) :: BspD, BspAug
  
  real :: x1, x2, y1, y2
  integer :: nx, ny
  integer :: N_subd
  integer :: i, j

  double precision :: gamma, hdiam, qexpo

  double precision, allocatable :: mysol(:)
  PetscScalar, pointer :: xx_v(:)

  Vec              :: sol          ! solution vector
  Vec              :: rhs
  Mat              :: Aps
  KSP              :: Solver
  PC               :: myPC
  PC               :: subPC
  KSP, allocatable :: subSolver(:)
  
  integer :: nlocal, first_local

  PetscErrorCode :: ierr
  PetscMPIInt    :: rank
  PetscBool      :: flg
  PetscErrorCode, parameter :: PETSC_OK = 0

  x1 = -1.0d0; x2 = 1.0d0; y1 = 1.0d0; y2 = 3.0d0
  nx = 10;   ny = 10;
  N_subd = 2

  call buildmesh2d(x1,x2,y1,y2,nx,ny,Th)
  call strip_dd(N_subd,Th,ddgeo)
  call plotmesh2d('box',Th)
  call plotinterface('interface',Th,ddgeo)

  hdiam = sqrt(Th%hx**2 + Th%hy**2)
  qexpo = 0.5d0
  gamma = 0.5d0 * (1.0d0 + hdiam**qexpo )
  !gamma = 0.5d0 * (1 + ((1.0d0/N_subd)*hdiam)**qexpo ) ! eta = 0
  !gamma = 0.5d0 * (1 + (hdiam/(1.0d0/N_subd))**qexpo ) ! eta = 1/h^2
  !gamma = 0.5d0 * (1 + (hdiam)**qexpo )                ! eta = 1/h

  !gamma = 1.0d0
  
  call fd_stiffness(Th,Asp)
  call build_res_opt(Th,ddgeo,res,resGAM)
  call build_loc_opt(Th,ddgeo,Asub,Asub2intf,Aintf)

  !! writing matrices into the file -- for two subdomain only
  call writespmat('scilab/A1',Asub(1))
  call writespmat('scilab/A2',Asub(2))
  !
  call writespmat('scilab/A1GAM',Asub2intf(1,1))
  call writespmat('scilab/A2GAM',Asub2intf(2,1))
  !
  call writespmat('scilab/AGAM',Aintf(1))
  !stop
  !! end writing matrices
  
  !  call showspmat(Asp)
  call sp_coo2csr(Asp,Bsp)
  ! zero based index for PETSC
  Bsp%i = Bsp%i - 1; Bsp%j = Bsp%j - 1

  call build_osm(ddgeo,osm)

  call build_augmatrix(osm, Asub, Asub2intf, Aintf, gamma, AaugD, AaugOff)

  do i=1, ddgeo%n_subd
     write(6,*) i, osm%subd_tot_dof(i)
  enddo

  call writespmat('AaugD',AaugD)
  call writespmat('AaugOff',AaugOff)

  call sp_coo2csr(AaugD, BspD)
  call sp_coo2csr(AaugD+AaugOff, BspAug)
  
  BspD%i = BspD%i - 1
  BspD%j = BspD%j - 1

  BspAug%i = BspAug%i - 1
  BspAug%j = BspAug%j - 1

  write(6,*) "START of PETSc"
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  !call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  call VecCreateSeq(PETSC_COMM_WORLD, BspAug%N_i ,rhs,ierr)
  call VecSet(rhs, 1.0d0,ierr)

  call VecCreateSeq(PETSC_COMM_WORLD, BspAug%N_i ,sol,ierr)

  call  MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, &
       BspAug%N_i, BspAug%N_j, BspAug%i, BspAug%j, BspAug%k, Aps, ierr)

  ! KSP solver
  call KSPcreate(PETSC_COMM_WORLD, Solver, ierr)
  call KSPSetOperators(Solver, Aps, Aps, ierr)
  !
  call KSPSetType(Solver, KSPRICHARDSON, ierr)
  ! Preconditioner
  call KSPGetPC(Solver, myPC, ierr)
  call PCSetType(myPC, PCBJACOBI, ierr)
  call PCBJacobiSetTotalBlocks(myPC, ddgeo%n_subd, osm%subd_tot_dof, ierr)
  ! End Precond.
  call KSPSetup(Solver, ierr)
  ! sub blocks precond
  allocate( SubSolver( ddgeo%n_subd ) )
  call PCBJacobiGetSubKSP(myPC, nlocal, first_local, SubSolver, ierr)

  do i=1, nlocal
     call KSPGetPC(SubSolver(i), SubPC, ierr)
     call PCSetType(SubPC, PCLU, ierr)
     ! call PCSetType(SubPC, PCICC, ierr)
     !call PCSetType(SubPC, PCCHOLESKY, ierr)
  enddo

  ! end sub precond
  !
     
  call KSPSetFromOptions(Solver, ierr)

  call KSPSolve(Solver, rhs, sol,ierr)

  call KSPView(Solver, PETSC_VIEWER_STDOUT_WORLD, ierr);
  !call VecView(sol, PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call VecView(rhs, PETSC_VIEWER_STDOUT_WORLD,ierr)

  !call MatAssemblyBegin(Aps,MAT_FINAL_ASSEMBLY,ierr)
  !call MatAssemblyEnd(Aps,MAT_FINAL_ASSEMBLY,ierr)
  !call MatView(Aps,PETSC_VIEWER_STDOUT_WORLD,ierr)

  allocate( mysol(BspAug%N_i) )
  call VecGetArrayF90(sol, xx_v, ierr)
  mysol = xx_v
  !call showvec( mysol )
  call VecRestoreArrayF90(sol,xx_v,ierr)

  call KSPGetIterationNumber(Solver, i)
  write(6,*) "size of system:", BspAug%N_i
  write(6,*) "hdiam:", hdiam, "gamma:", gamma
  write(6,*) "number of iterations: ", i
  write(6,*) "hdim, hdim**qexp: ", hdiam, hdiam**qexpo
  call KSPDestroy(Solver,ierr)
  call VecDestroy(sol,ierr)
  call VecDestroy(rhs,ierr)
  call MatDestroy(Aps,ierr)

  call PetscFinalize(ierr) 

contains

  !! -------------------------------
  !> Displaying a Matrix in the shell
  !! @param [in] M matrix of size i,j
  subroutine ShowMat(M)
    real, intent(in), dimension(:,:) :: M
    integer :: dimx, dimy, i
    character(len=80) :: lb

    dimx = size(M, dim=1)
    dimy = size(M, dim=2)

    do i = 1, dimx
       lb = "('[',I2,']', 1000F8.2)"
       write(*,lb) i, M(i,:)
    enddo
    write(*,"(/)")

  end subroutine ShowMat

  subroutine ShowVec(V)
    double precision, intent(in), dimension(:) :: V
    integer :: dimx, i
    character(len=80) :: lb
    dimx = size(V, dim=1)

    lb = "('[',I4,']', 1000F8.4)"
    do i=1, dimx
       write(*,lb) i, V(i)
    enddo
    write(*,"(/)")
  end subroutine ShowVec

end program main_petsc


