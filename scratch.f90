  !> scratch.f90

program scratch

  use fd

  implicit none

  type(mesh2d) :: Th
  type(dd1d)   :: ddgeo
  type(spmat)  :: Asp
  type(spmat), allocatable  :: res(:), resGAM(:)
  type(spmat), allocatable  :: Asub(:), Asub2intf(:,:), Aintf(:)
  type(spmat) :: Csp, Dsp

  double precision, allocatable :: vecX(:)
  
  real :: x1, x2, y1, y2
  integer :: nx, ny
  integer :: N_subd
  integer :: i, j

  x1 = -1.0; x2 = 1.0; y1 = 1.0; y2 = 3.0
  nx = 20*16;   ny = 20*16
  N_subd = 9
  !
  ! call buildspmat(Csp,100,120,50)
  ! do i=1, 3
  !    call insert2spmat(Csp,i,i+1,0.33*i)
  !    call insert2spmat(Csp,i,i+2,0.5*i)
  !    call insert2spmat(Csp,i,i+1,0.33*i)
  ! end do

  ! call copysp2sp(Csp,Dsp)
  ! call copysp2sp(Csp,Dsp)
  ! call copysp2sp(Csp,Dsp)

  ! call showspmat(Csp)
  ! call sum_redundant(Dsp,Dsp)
  ! call showspmat(Dsp)
  ! call sum_redundant(Dsp,Dsp)
  ! call showspmat(Dsp)

  ! stop

  !
  call buildmesh2d(x1,x2,y1,y2,nx,ny,Th)
  call strip_dd(N_subd,Th,ddgeo)
  call plotmesh2d('box',Th)
  call plotinterface('interface',Th,ddgeo)

  call fd_stiffness(Th,Asp)
  call build_res_opt(Th,ddgeo,res,resGAM)
  call build_loc_opt(Th,ddgeo,Asub,Asub2intf,Aintf)

  ! allocate( vecX( Asp%N_j ) )
  ! vecX = 1.0d0
  ! vecX = Asp * VecX
  ! call showvec( vecX )
  ! stop
  
  ! write matrices for vertification
  call writespmat('scilab/A1',Asub(1))
  call writespmat('scilab/A1GAM',Asub2intf(1,1))
  call writespmat('scilab/AGAM',Aintf(1))

  call writespmat('scilab/A2',Asub(2))
  call writespmat('scilab/A2GAM',Asub2intf(2,1))
  stop
  ! write matrices for vertification
  call writespmat('scilab/A2',Asub(2))
  call writespmat('scilab/A2GAM1',Asub2intf(2,1))
  call writespmat('scilab/A2GAM2',Asub2intf(2,2))
  call writespmat('scilab/AGAM1',Aintf(1))
  call writespmat('scilab/AGAM2',Aintf(2))
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
end program scratch

