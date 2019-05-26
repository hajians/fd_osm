!> test1 program checks whether or not the 9-point stencil restriction
!! generates a coarse stiffness matrix on a coarser grid.
program test1

  use fd
  use coarsetools
  
  implicit none

  type(mesh2d)  :: Th, ThC     ! Th for fine mesh, ThC for coarse mesh
  type(dd1d)    :: ddgeo
  type(coarseT) :: mycoarse    ! coarse tools
  type(spmat)   :: Asp, AspC   ! Asp for fine mesh, AspC for coarse mesh

  real    :: x1, x2, y1, y2        ! coordinates of the box
  integer :: nx, ny, nxC, nyC      ! number of nodes in fine and coarse meshes
  integer :: FCfactorX, FCfactorY  ! fine/coarse factor
  integer :: N_subd                ! number of subdomains

  x1 = -1.0; x2 = 1.0; y1 = 1.0; y2 = 3.0

  FCfactorX = 2; FCfactorY = 2
  
  nx = 12; ny = 12; nxC = nx / FCfactorX; nyC = ny / FCfactorY;

  N_subd = 2

  call buildmesh2d(x1,x2,y1,y2,nx,ny,Th)
  call buildmesh2d(x1,x2,y1,y2,nxC,nyC,ThC)
  call strip_dd(N_subd,Th,ddgeo)
  call buildcoarsetools(ThC,Th,mycoarse)
  
  call plotmesh2d('fine',Th)
  call plotmesh2d('coarse',ThC)

  call fd_stiffness(Th,Asp)
  call fd_stiffness(ThC,AspC)

  write(6,*) "%% writing stiffness matrices into file"
  call writespmat('scilab/fineStiffness',Asp)
  call writespmat('scilab/coarseStiffness',AspC)

  write(6,*) "%% writing prolongation matrix"
  call writespmat('scilab/prolongation', mycoarse%prolongation)

  call ShowVec( real(mycoarse%Cv2Fv) )

  call showspmat( mycoarse%prolongation )


contains

  !> Displaying a Matrix in the shell
  !! @param [in] M matrix of size i,j
  subroutine ShowMat(M)
    double precision, intent(in), dimension(:,:) :: M
    integer :: dimx, dimy, i
    character(len=80) :: lb

    dimx = size(M, dim=1)
    dimy = size(M, dim=2)

    do i = 1, dimx
       lb = "('[',I3,']', 1000F8.2)"
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

end program test1
