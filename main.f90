!> main.f90 This is the main file of this program.  It generates
!> matrices necessary for domain decomposition purposes
program main

  use fd

  implicit none

  type(mesh2d) :: Th
  type(dd1d)   :: ddgeo
  type(spmat)  :: Asp
  type(spmat), allocatable  :: res(:), resGAM(:)
  type(spmat), allocatable  :: Asub(:), Asub2intf(:,:), Aintf(:)
  ! type(spmat) :: Csp, Dsp

  ! double precision, allocatable :: vecX(:)
  
  real :: x1, x2, y1, y2
  integer :: nx, ny
  integer :: N_subd
  ! integer :: i, j

  x1 = -1.0; x2 = 1.0; y1 = 1.0; y2 = 3.0
  nx = 20*1;   ny = 20*1
  N_subd = 2

  write(*,*) "nx: ", nx, " ny: ", ny
  
  call buildmesh2d(x1,x2,y1,y2,nx,ny,Th)
  call strip_dd(N_subd,Th,ddgeo)
  call plotmesh2d('box',Th)
  call plotinterface('interface',Th,ddgeo)

  call fd_stiffness(Th,Asp)
  call build_res_opt(Th,ddgeo,res,resGAM)
  call build_loc_opt(Th,ddgeo,Asub,Asub2intf,Aintf)

  write (*,*) 'saving matrices in scilab/data/ folder'
  ! write matrices for vertification
  call writespmat('scilab/data/A1',Asub(1))
  call writespmat('scilab/data/A1GAM',Asub2intf(1,1))
  call writespmat('scilab/data/AGAM',Aintf(1))

  call writespmat('scilab/data/A2',Asub(2))
  call writespmat('scilab/data/A2GAM',Asub2intf(2,1))
  stop
  ! write matrices for vertification
  call writespmat('scilab/data/A2',Asub(2))
  call writespmat('scilab/data/A2GAM1',Asub2intf(2,1))
  call writespmat('scilab/data/A2GAM2',Asub2intf(2,2))
  call writespmat('scilab/data/AGAM1',Aintf(1))
  call writespmat('scilab/data/AGAM2',Aintf(2))
contains

  !! -------------------------------
  !> Displaying a Matrix in the shell
  !! @param [in] M matrix of size i,j
  subroutine ShowMat(M)
    double precision, intent(in), dimension(:,:) :: M
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
  !! -------------------------------
end program main


