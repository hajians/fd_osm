!> @file mesh2d.f90 produce an structured mesh over a box
!!

module geo2d
  implicit none

  !> mesh2d is an object that represents a 2d mesh over a box
  !! in 2D. It is useful for a finite difference method.
  type mesh2d

     !> initialization flag
     logical :: init = .false.

     !> dimension of the domain. It will be set to 2.
     integer :: n_dim

     !> total number of vertices
     integer :: n_vt, n_vtx, n_vty

     !> number of elements
     integer :: nx, ny

     !> mesh size
     real    :: hx, hy

     !> vertices to coordinate map
     real,    allocatable, dimension(:,:) :: v2corr

     !> vertices to neighboring vertices map
     integer, allocatable, dimension(:,:) :: v2v

     !> coordinates of the domain
     real,    allocatable, dimension(:,:) :: omega

  end type mesh2d

  !> vertex represent a subdomain or interface
  !!  by storing the corresponding DOF's index.
  type vertex
     integer :: n_index
     integer, allocatable, dimension(:) :: index
  end type vertex

  !> dd1d is an object which represent a 1d domain decomposition,
  !! in other words we decompose the domain in one direction only;
  !! decompostion of the domain into strips.
  type dd1d

     !> initialization flag
     logical :: init = .false.

     !> number of subdomains
     integer :: n_subd

     !> number of interfaces
     integer :: n_intf

     !> sets of subdomain vertices
     type(vertex), allocatable, dimension(:) :: subdm

     !> sets of interfaces
     type(vertex), allocatable, dimension(:) :: interf

     !> subdomain to interface map
     integer, allocatable, dimension(:,:)    :: subd2intf

     !> interface to subdomain map
     integer, allocatable, dimension(:,:)    :: intf2subd

     !> global vertex index to interface map. If a vertex is not on the
     !! interface we assign 0.
     integer*2, allocatable, dimension(:)    :: v2intf

     !> global vertex index to subdomain and local map.
     !! If a vertex does not belong to a subdomain, we assign 0.
     integer, allocatable, dimension(:,:)    :: glob2loc
     
  end type dd1d


contains

  !> builds the 2d mesh, Th, using the coordinates of the
  !! domain: [x1,x2]*[y1,y2] and with # mesh elements
  !! nx and ny. It is supposed x2>x1, y2>y1.
  !! We do NOT save the boundary vertices.
  subroutine buildmesh2d(x1,x2,y1,y2,nx,ny,Th)
    real,         intent(in)    :: x1, x2, y1, y2
    integer,      intent(in)    :: nx, ny
    type(mesh2d), intent(inout) :: Th

    integer :: i, j, k

    if ( Th%init.eqv..true. ) then
       write(6,*) "error in buildmesh2d: Th%init"
    endif

    Th%n_dim = 2                ! dimension is 2

    Th%nx = nx
    Th%ny = ny

    Th%hx = (x2-x1)/nx          ! mesh size - x
    Th%hy = (y2-y1)/ny          ! mesh size - y

    Th%n_vtx = nx - 1           ! # of vertices in x
    Th%n_vty = ny - 1           ! # of vertices in y
    Th%n_vt  = Th%n_vtx * Th%n_vty ! total # of vertices

    allocate( Th%omega(2,2) )
    Th%omega(1,:) = (/ x1, x2 /)
    Th%omega(2,:) = (/ y1, y2 /)

    allocate( Th%v2corr(Th%n_vt, Th%n_dim) )

    allocate( Th%v2v(Th%n_vt, 4) ) ! FOUR neighbors

    do j=1, Th%n_vty
       do i=1, Th%n_vtx
          k = i + (j-1)*Th%n_vtx ! vertex index
          Th%v2corr(k,:) = (/ x1 + i*Th%hx , x2 + j*Th%hy /)
       enddo
    enddo

    !! building the adjacency map for vertices
    Th%v2v = -1

    do j=1, Th%n_vty            ! We start with "core" nodes.
       do i=1, Th%n_vtx         ! Same here.
          k = i + (j-1)*Th%n_vtx  ! vertex index

          if ( i.le.(Th%n_vtx-1) ) then
             Th%v2v(k,2) = k + 1
          endif
          if ( i.ge.2 ) then
             Th%v2v(k,1) = k - 1
          endif
          if( j.le.(Th%n_vty-1) ) then
             Th%v2v(k,4) = k + Th%n_vtx
          endif
          if ( j.ge.2 ) then
             Th%v2v(k,3) = k - Th%n_vtx
          endif

       end do
    end do

    Th%init = .true.
  end subroutine buildmesh2d

  !> plot a mesh2d object, Th, in a file.
  subroutine plotmesh2d(filename, Th)
    character(len=*),  intent(in) :: filename
    type(mesh2d),      intent(in) :: Th

    integer :: i
    integer :: myunit = 10

    open(myunit, file=filename//'.gnu')
    do i=1, Th%n_vt
       write(myunit,*) Th%v2corr(i,:)
    enddo
    close(myunit)
  end subroutine plotmesh2d

  !> strip decomposition of the domain.
  !! it takes # of subdomain, N_subd, and mesh, Th.
  !! it builds the domain decomposition maps, geo_dd 
  subroutine strip_dd(N_subd,Th,geo_dd)
    integer,      intent(inout) :: N_subd
    type(mesh2d), intent(in)    :: Th
    type(dd1d),   intent(inout) :: geo_dd

    integer :: max_strip        ! maximum number of strips
    integer :: dummy
    integer :: subd_strip
    integer :: subd_extra
    integer :: subd, intf, vtx
    integer :: counter
    integer :: i, j, k

    if ( Th%init.eqv..false. .or. geo_dd%init.eqv..true. ) then
       write(6,*) "error in strip_dd: "
       stop
    endif

    !! structure of the domain decomposition
    !! | Omega_1 | Omega_2 | ................. | Omega_(N_subd) |

    max_strip = ceiling( real(Th%n_vtx - 2)/2 )

    if ( N_subd > (max_strip+1) ) then
       N_subd = max_strip+1
       write(6,*) "warning in strip_dd: N_subd is modified to", N_subd
    endif

    geo_dd%n_subd = N_subd
    geo_dd%n_intf = N_subd - 1

    allocate( geo_dd%subdm( N_subd )   ) ! allocate subdomains
    allocate( geo_dd%interf(N_subd-1 ) ) ! allocate interfaces
    allocate( geo_dd%v2intf(Th%n_vt) )   ! allocate v2intf
    allocate( geo_dd%glob2loc(Th%n_vt, N_subd) )

    geo_dd%v2intf   = 0         ! initialize v2intf
    geo_dd%glob2loc = 0         ! initialize glob2loc

    dummy = Th%n_vtx - geo_dd%n_intf 

    subd_strip = dummy/N_subd        ! 
    subd_extra = mod(dummy,N_subd)   !

    ! first subdomain
    geo_dd%subdm(1)%n_index = (subd_strip+subd_extra)*Th%n_vty
    allocate( geo_dd%subdm(1)%index( geo_dd%subdm(1)%n_index) )

    counter = 1                 !
    !
    do j=1, Th%n_vty
       do i=1, (subd_strip+subd_extra)
          k = i + (j-1)*Th%n_vtx
          geo_dd%subdm(1)%index( counter ) = k
          counter = counter + 1
       end do
    enddo
    vtx = 1 + (subd_strip+subd_extra)
    ! end of subdomain one
    !
    !
    ! loop over subdomains > 1
    do subd = 2, N_subd
       geo_dd%subdm(subd)%n_index = subd_strip*Th%n_vty
       allocate( geo_dd%subdm(subd)%index( geo_dd%subdm(subd)%n_index) )
       counter = 1                 !

       do j=1, Th%n_vty
          do i=1, subd_strip
             k = (i+vtx) + (j-1)*Th%n_vtx
             geo_dd%subdm(subd)%index( counter ) = k
             counter = counter + 1
          end do
       enddo
       vtx = vtx + 1 + subd_strip
    enddo


    ! loop for interfaces
    vtx = 1 + (subd_strip+subd_extra)
    do intf=1, geo_dd%n_intf
       geo_dd%interf(intf)%n_index = Th%n_vty
       allocate( geo_dd%interf(intf)%index( geo_dd%interf(intf)%n_index) )
       counter = 1
       do j=1, Th%n_vty
          k = vtx + (j-1)*Th%n_vtx
          geo_dd%interf(intf)%index(counter) = k
          counter = counter + 1
       enddo
       vtx = vtx + subd_strip + 1
    enddo
    ! end loop for interfaces

    ! subdomain-to-interface map and interface-to-subdomain map
    allocate( geo_dd%subd2intf( N_subd, 2 ) )      ! maximum two adjacent interface
    allocate( geo_dd%intf2subd( geo_dd%n_intf, 2 ))! two adjacent subdomain

    geo_dd%subd2intf      = -1
    !
    geo_dd%subd2intf(1,2) = 1
    geo_dd%subd2intf(N_subd,1) = geo_dd%n_intf
    !
    do i=2, N_subd-1
       geo_dd%subd2intf(i,1) = i-1
       geo_dd%subd2intf(i,2) = i
    enddo

    do i=1, geo_dd%n_intf
       geo_dd%intf2subd(i,1) = i
       geo_dd%intf2subd(i,2) = i+1
    enddo

    do intf=1, geo_dd%n_intf
       do j=1, geo_dd%interf(intf)%n_index
          k = geo_dd%interf(intf)%index(j) ! k is the vertex index on interface intf
          geo_dd%v2intf(k) = intf
       enddo
    enddo

    do subd=1, N_subd
       do i=1, geo_dd%subdm(subd)%n_index ! local  index
          j = geo_dd%subdm(subd)%index(i) ! global index
          geo_dd%glob2loc(j,subd) = i     ! fill the map
       enddo
    enddo

    geo_dd%init = .true.
  end subroutine strip_dd

  !> plot interface.
  subroutine plotinterface(filename, Th, geo_dd)
    character(len=*),  intent(in) :: filename
    type(mesh2d),      intent(in) :: Th
    type(dd1d),        intent(in) :: geo_dd

    integer :: i, j, vertex
    integer :: myunit = 10

    open(myunit, file=filename//'.gnu')
    do i=1, geo_dd%n_intf
       do j=1, geo_dd%interf(i)%n_index
          vertex = geo_dd%interf(i)%index(j)
          write(myunit,*) Th%v2corr(vertex,:)
       enddo
    enddo
    close(myunit)
  end subroutine plotinterface

end module geo2d
