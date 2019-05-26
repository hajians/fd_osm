!> fd is a module that contains subroutines for
!! building finite difference discretizations and
!! domain decomposition components.
module fd
  use geo2d
  use sparsemat

  implicit none

  !> OSMopt contains operators to build the augmented matrix
  !! associated to OSM on a STRIP.
  type OSMopt
     
     !! number of subdomain
     integer :: N_subd
     
     !! total number of primal unknowns
     integer :: primal_tot_dofs

     !! total number of interface unknowns
     integer :: interf_tot_dofs
     
     !! total number of unknowns
     integer :: tot_dofs

     !! total number of unknowns per subdomain
     integer, allocatable :: subd_tot_dof(:)
     
     !! subdomain primal unknowns to global dofs
     type(vertex), allocatable :: subd_primal2Gdof(:)

     !! subdomain interface unknowns to global dofs.
     !! Except first and last subdomain, other subdomains
     !! have two interfaces.
     type(vertex), allocatable :: subd_interf2Gdof(:,:)

     !! all subdomain global dofs to mesh vertices
     integer, allocatable :: Gdof2vert(:)
  end type OSMopt
  
contains

  !> fd_stiffness builds stiffness matrix
  !! based on the triangulation Th and writes
  !! in Asp, a sparse matrix.
  !! Asp should NOT had been initialized.
  subroutine fd_stiffness(Th,Asp)
    type(mesh2d), intent(in)    :: Th
    type(spmat),  intent(inout) :: Asp

    integer :: n_dof            ! total number of dofs

    double precision    :: hx               ! hx is the mesh size in x direction
    double precision    :: hy               ! hy is the mesh size in y direction

    ! in order to save time we compute -1/hx^2, +2/hx^2, -1/hy^2, +2/hy^2
    double precision    :: mHX, p2HX, mHY, p2HY


    
    integer :: vert             ! vertex index
    integer :: adjVert          ! adjacent vertex index

    integer :: i

    if ( (Asp%init.eqv..true.) .or. (Th%init.eqv..false.) ) then
       write(6,*) "error in fd_stiffness", Th%init, Asp%init
       stop
    endif

    n_dof = Th%n_vt             ! number of dofs
    hx    = Th%hx               ! hx
    hy    = Th%hy               ! hy    

    mHX  = -1.0 / hx**2 ; mHY  = -1.0 / hy**2
    p2HX = +2.0 / hx**2 ; p2HY = +2.0 / hy**2


    
    call buildspmat(Asp,n_dof,n_dof, 5*n_dof) ! # non-zero entries ~ 5*n_dof

    do vert=1, Th%n_vt
       call insert2spmat(Asp,vert,vert, p2HX + p2HY)
       ! adjacent filling
       do i=1, size(Th%v2v, 2)
          adjVert = Th%v2v(vert, i)
          ! if adjVert is not on the boundary
          if (adjVert.ne.-1) then
             !! THIS IS NOT NICE: WE NEED TO RELY ON NUMBERING OF Th.
             if ( (i.le.2) ) then
                call insert2spmat(Asp,vert,adjVert, mHX)             
             else
                call insert2spmat(Asp,vert,adjVert, mHY)             
             endif
             !! END OF NOT-NICENESS.
          endif
          ! end adding adjacent
       end do
    enddo

    ! clean up zero entries
    call deleteblankspmat(Asp)

  end subroutine fd_stiffness

  !> build_res_opt builds the restriction operators
  !! and write them in the sparse matrices.
  !! res is an array of sparse matrices for subdomains
  !! and resGam is an array of sparse matrices for interfaces.
  subroutine build_res_opt(Th, geo_dd,res,resGam)
    type(mesh2d), intent(in)                 :: Th
    type(dd1d),   intent(in)                 :: geo_dd
    type(spmat),  intent(inout), allocatable :: res(:)
    type(spmat),  intent(inout), allocatable :: resGAM(:)

    integer :: N_subd
    integer :: N_intf
    integer :: i, j

    if ( (geo_dd%init.eqv..false.) .or. allocated(res) .or. allocated(resGAM) ) then
       write(6,*) "error in build_res_opt"
       stop
    end if

    N_subd = geo_dd%n_subd
    N_intf = geo_dd%n_intf
    allocate( res( N_subd) )
    allocate( resGAM( N_intf) )

    !! building restriction for the subdomains
    do i=1, N_subd
       call buildspmat(res(i), geo_dd%subdm(i)%n_index, Th%n_vt, geo_dd%subdm(i)%n_index)
       do j=1, geo_dd%subdm(i)%n_index
          call insert2spmat(res(i), j, geo_dd%subdm(i)%index(j), 1.0d0)
       end do
       call deleteblankspmat(res(i))
    end do

    !! building restriction for the interfaces
    do i=1, N_intf
       call buildspmat(resGAM(i), geo_dd%interf(i)%n_index, Th%n_vt, &
            geo_dd%interf(i)%n_index )
       do j=1, geo_dd%interf(i)%n_index
          call insert2spmat(resGAM(i), j, geo_dd%interf(i)%index(j), 1.0d0)
       enddo
       call deleteblankspmat(resGAM(i))       
    end do

  end subroutine build_res_opt

  !> build_loc_opt builds local stiffness matrices
  !! and matrices between subdomains and interfaces.
  subroutine build_loc_opt(Th,geo_dd,Asub,Asub2intf,Aintf)
    type(mesh2d), intent(in)                 :: Th
    type(dd1d),   intent(in)                 :: geo_dd
    type(spmat),  intent(inout), allocatable :: Asub(:)
    type(spmat),  intent(inout), allocatable :: Asub2intf(:,:)
    type(spmat),  intent(inout), allocatable :: Aintf(:)

    integer :: N_subd
    integer :: N_intf
    integer :: sub_dof
    integer :: intf_dof
    integer :: subd
    integer :: intf
    integer :: vert, adjVert
    integer :: loc_dof, adj_loc_dof

    integer :: i, j

    double precision    :: hx, hy
    double precision    :: mHX, p2HX, mHY, p2HY
    double precision    :: Eta
    
    if ( (geo_dd%init.eqv..false.) .or. &
         allocated(Asub) .or. allocated(Aintf) .or. allocated(Asub2intf)  ) then
       write(6,*) "error in build_loc_opt"
       stop
    end if

    N_subd = geo_dd%n_subd
    N_intf = geo_dd%n_intf

    allocate( Asub(N_subd) )
    allocate( Aintf(N_intf) )
    allocate( Asub2intf(N_subd, N_intf) )

    hx    = Th%hx               ! hx
    hy    = Th%hy               ! hy    

    mHX  = -1.0 / hx**2 ; mHY  = -1.0 / hy**2
    p2HX = +2.0 / hx**2 ; p2HY = +2.0 / hy**2
    !Eta = 1.0 / sqrt(hx**2 + hy**2)
    Eta = 0.0
    
    !! start building the local solves
    do subd=1, N_subd
       sub_dof = geo_dd%subdm(subd)%n_index
       call buildspmat(Asub(subd), sub_dof, sub_dof, 5*sub_dof)
       !
       do i=1, sub_dof
          vert    = geo_dd%subdm(subd)%index(i)
          loc_dof = geo_dd%glob2loc(vert,subd) ! this is redundant
          !
          call insert2spmat(Asub(subd),loc_dof,loc_dof, Eta + p2HX + p2HY)
          ! Adjacent nodes - filling
          do j=1, size(Th%v2v,2)               ! max. four neighbor per vertex
             adjVert = Th%v2v(vert,j)
             !
             if ( adjVert.ne.-1 ) then         ! if not on the boundary of domain
                adj_loc_dof = geo_dd%glob2loc(adjVert,subd)
                if ( adj_loc_dof.ne.0 ) then   ! if not on the interface
                   !
                   if (j.le.2) then            ! if neighbor is on the x-axis
                      !! THIS IS NOT NICE: WE NEED TO RELY ON NUMBERING OF Th.
                      call insert2spmat(Asub(subd),loc_dof,adj_loc_dof,mHX)
                   else
                      call insert2spmat(Asub(subd),loc_dof,adj_loc_dof,mHY)
                   endif
                   !
                else

                endif
             endif
             !
          enddo
          !
       enddo
       !
       call deleteblankspmat(Asub(subd))
    enddo                       ! end of subdomain loop

    !! building subdomain-to-interface matrices
    do intf=1, N_intf
       intf_dof = geo_dd%interf(intf)%n_index ! # of dofs on the interface, intf
       do j=1, size(geo_dd%intf2subd,2)       ! loop over neighboring subdomains
          subd = geo_dd%intf2subd(intf,j)     ! neighboring subdomain
          sub_dof = geo_dd%subdm(subd)%n_index! number of unknowns in subdomain subd
          !! build matrix
          call buildspmat(Asub2intf(subd,intf), sub_dof, intf_dof, 5*sub_dof)
          !!
          do i=1, intf_dof
             vert     = geo_dd%interf(intf)%index(i)

             adjVert     = Th%v2v(vert,j)                ! j=1 means left, j=2 means right
             adj_loc_dof = geo_dd%glob2loc(adjVert,subd)
             if ( adj_loc_dof.eq.0 ) then
                write(6,*) "bug in build_loc_opt: subdomain-to-interface"
                stop
             endif
             ! insert
             call insert2spmat(Asub2intf(subd,intf), adj_loc_dof, i, mHX)
             ! end of insertion
          enddo
          !!
          call deleteblankspmat(Asub2intf(subd,intf))
       enddo                                  ! end of visiting subdomain

    enddo                       ! end of visiting the interfaces

    !! filling the interface matrices
    do intf=1, N_intf
       intf_dof = geo_dd%interf(intf)%n_index ! # of dofs on the interface, intf
       call buildspmat(Aintf(intf), intf_dof, intf_dof, 5*intf_dof)
       !!
       do i=1, intf_dof
          call insert2spmat(Aintf(intf),i,i, Eta + p2HX + p2HY)
          if (i>1) then
             call insert2spmat(Aintf(intf),i,i-1, mHY)
          endif
          if (i<intf_dof) then
             call insert2spmat(Aintf(intf),i,i+1, mHY)
          endif
       end do
       !!
       call deleteblankspmat(Aintf(intf))
    enddo

  end subroutine build_loc_opt

  !> build osm_opt from geo_dd structure for a domain decomposition
  !! in a strip.
  subroutine build_osm(geo_dd,osm_opt)
    type(dd1d),   intent(in)         :: geo_dd
    type(OSMopt), intent(inout)      :: osm_opt

    integer :: i, j, k
    integer :: dummy
    integer :: tot_dofs

    tot_dofs = 0
    !
    do i=1, geo_dd%n_subd
       tot_dofs = tot_dofs + geo_dd%subdm(i)%n_index
    enddo

    osm_opt%primal_tot_dofs = tot_dofs

    tot_dofs = 0
    !
    do i=1, geo_dd%n_intf
       tot_dofs = tot_dofs + geo_dd%interf(i)%n_index
    enddo
    
    osm_opt%interf_tot_dofs = 2 * tot_dofs

    osm_opt%tot_dofs = osm_opt%interf_tot_dofs + osm_opt%primal_tot_dofs

    osm_opt%N_subd = geo_dd%n_subd
    
    allocate( osm_opt%subd_primal2Gdof( geo_dd%n_subd ) )
    allocate( osm_opt%subd_interf2Gdof( geo_dd%n_subd, 2 ) ) ! two
                                                             ! interface
                                                             ! per
                                                             ! subdomain.
    allocate( osm_opt%Gdof2vert( osm_opt%tot_dofs ) )
    osm_opt%Gdof2vert = 0

    allocate( osm_opt%subd_tot_dof( geo_dd%n_subd ) )


    

    !! filling gdofs
    dummy = 0
    i = 1                       ! subdomain one
    !
    osm_opt%subd_primal2Gdof(i)%n_index = geo_dd%subdm(i)%n_index
    allocate( osm_opt%subd_primal2Gdof(i)%index( geo_dd%subdm(i)%n_index ) )
    osm_opt%subd_primal2Gdof(i)%index = &
         (/ (k, k=1, geo_dd%subdm(i)%n_index) /)  + dummy
    dummy = osm_opt%subd_primal2Gdof(i)%index( geo_dd%subdm(i)%n_index )
    osm_opt%Gdof2vert( osm_opt%subd_primal2Gdof(i)%index ) = &
         geo_dd%subdm(i)%index
    !
    j = 2
    osm_opt%subd_interf2Gdof(i,j)%n_index = &
         geo_dd%interf(i-1+j-1)%n_index
    allocate( osm_opt%subd_interf2Gdof(i,j)%index( &
         geo_dd%interf(i-1+j-1)%n_index ) )
    osm_opt%subd_interf2Gdof(i,j)%index = &
         (/ (k, k=1, geo_dd%interf(i-1+j-1)%n_index) /) + dummy
    dummy = osm_opt%subd_interf2Gdof(i,j)%index( geo_dd%interf(i-1+j-1)%n_index )
    !
    osm_opt%Gdof2vert( osm_opt%subd_interf2Gdof(i,j)%index ) = &
         geo_dd%interf(i-1+j-1)%index    
    !
    ! other subdomains
    do i=2, (geo_dd%n_subd-1)
       osm_opt%subd_primal2Gdof(i)%n_index = geo_dd%subdm(i)%n_index
       allocate( osm_opt%subd_primal2Gdof(i)%index( geo_dd%subdm(i)%n_index ) )
       osm_opt%subd_primal2Gdof(i)%index = &
            (/ (k, k=1, geo_dd%subdm(i)%n_index) /)  + dummy
       dummy = osm_opt%subd_primal2Gdof(i)%index( geo_dd%subdm(i)%n_index )
       !
       osm_opt%Gdof2vert( osm_opt%subd_primal2Gdof(i)%index ) = &
            geo_dd%subdm(i)%index
       !
       do j=1, 2                ! two faces
          osm_opt%subd_interf2Gdof(i,j)%n_index = &
               geo_dd%interf(i-1+j-1)%n_index
          !
          allocate( osm_opt%subd_interf2Gdof(i,j)%index( &
               geo_dd%interf(i-1+j-1)%n_index ) )
          !
          osm_opt%subd_interf2Gdof(i,j)%index = &
          (/ (k, k=1, geo_dd%interf(i-1+j-1)%n_index) /) + dummy
          dummy = osm_opt%subd_interf2Gdof(i,j)%index( geo_dd%interf(i-1+j-1)%n_index )
          !
          osm_opt%Gdof2vert( osm_opt%subd_interf2Gdof(i,j)%index ) = &
               geo_dd%interf(i-1+j-1)%index    
       enddo
    enddo
    ! last subdomain
    i = geo_dd%n_subd                       ! subdomain one
    !
    osm_opt%subd_primal2Gdof(i)%n_index = geo_dd%subdm(i)%n_index
    allocate( osm_opt%subd_primal2Gdof(i)%index( geo_dd%subdm(i)%n_index ) )
    osm_opt%subd_primal2Gdof(i)%index = &
         (/ (k, k=1, geo_dd%subdm(i)%n_index) /)  + dummy
    dummy = osm_opt%subd_primal2Gdof(i)%index( geo_dd%subdm(i)%n_index )
    osm_opt%Gdof2vert( osm_opt%subd_primal2Gdof(i)%index ) = &
         geo_dd%subdm(i)%index
    !
    j = 1
    osm_opt%subd_interf2Gdof(i,j)%n_index = &
         geo_dd%interf(i-1+j-1)%n_index
    allocate( osm_opt%subd_interf2Gdof(i,j)%index( &
         geo_dd%interf(i-1+j-1)%n_index ) )
    osm_opt%subd_interf2Gdof(i,j)%index = &
         (/ (k, k=1, geo_dd%interf(i-1+j-1)%n_index) /) + dummy
    dummy = osm_opt%subd_interf2Gdof(i,j)%index( geo_dd%interf(i-1+j-1)%n_index )
    !
    osm_opt%Gdof2vert( osm_opt%subd_interf2Gdof(i,j)%index ) = &
         geo_dd%interf(i-1+j-1)%index    

    do i=1, geo_dd%n_subd
       osm_opt%subd_tot_dof(i) = osm_opt%subd_primal2Gdof(i)%n_index
       if (i>1) then
          osm_opt%subd_tot_dof(i) = osm_opt%subd_tot_dof(i) + &
               osm_opt%subd_interf2Gdof(i,1)%n_index
       endif
       if (i < geo_dd%n_subd) then
          osm_opt%subd_tot_dof(i) = osm_opt%subd_tot_dof(i) + &
               osm_opt%subd_interf2Gdof(i,2)%n_index
       endif
    enddo

  end subroutine build_osm

  !> using osm_opt and block matrices Asub, Asub2intf, Aintf builds
  !! augmented matrix Aaug
  subroutine build_augmatrix(osm_opt, Asub, Asub2intf, Aintf, gamma, AaugD, AaugOff)
    type(OSMopt), intent(in)    :: osm_opt
    type(spmat),  intent(in)    :: Asub(:)
    type(spmat),  intent(in)    :: Asub2intf(:,:)
    type(spmat),  intent(in)    :: Aintf(:)
    type(spmat),  intent(inout) :: AaugD, AaugOff

    double precision, intent(in) :: gamma
    
    integer :: i

    if ( (AaugD%init.eqv..true.) .or. (AaugOff%init.eqv..true.) ) then
       write(6,*) "error in build_augmatrix: "
       stop
    endif

    call buildspmat(AaugD, osm_opt%tot_dofs, osm_opt%tot_dofs, 5*osm_opt%tot_dofs)
    call buildspmat(AaugOff, osm_opt%tot_dofs, osm_opt%tot_dofs, 5*osm_opt%tot_dofs)

    !! fill the SUBDOMAIN diagonals
    do i=1, osm_opt%N_subd
       call insertblock2spmat(AaugD, &
            osm_opt%subd_primal2Gdof(i)%index(1), osm_opt%subd_primal2Gdof(i)%index(1), &
            Asub(i) )
       !
       if (i > 1) then
          !! here we should define $gamma$
          call insertblock2spmat(AaugD, &
               osm_opt%subd_interf2Gdof(i,1)%index(1), &
               osm_opt%subd_interf2Gdof(i,1)%index(1), &
               gamma * Aintf(i-1) )     ! 1 stands for left
          !! sub2interf
          call insertblock2spmat(AaugD, &
               osm_opt%subd_primal2Gdof(i)%index(1), &
               osm_opt%subd_interf2Gdof(i,1)%index(1), &
               Asub2intf(i,i-1) ) ! 1 stands for left
          call insertblock2spmat(AaugD, &
               osm_opt%subd_interf2Gdof(i,1)%index(1), &
               osm_opt%subd_primal2Gdof(i)%index(1), &
               SPtranspose(Asub2intf(i,i-1)) ) ! 1 stands for left
       endif
       !! here we should define $gamma$
       if (i < osm_opt%N_subd) then
          call insertblock2spmat(AaugD, &
               osm_opt%subd_interf2Gdof(i,2)%index(1), &
               osm_opt%subd_interf2Gdof(i,2)%index(1), &
               gamma * Aintf(i) )     ! 2 stands for right
          call insertblock2spmat(AaugD, &
               osm_opt%subd_primal2Gdof(i)%index(1), &
               osm_opt%subd_interf2Gdof(i,2)%index(1), &
               Asub2intf(i,i) ) ! 2 stands for right
          call insertblock2spmat(AaugD, &
               osm_opt%subd_interf2Gdof(i,2)%index(1), &
               osm_opt%subd_primal2Gdof(i)%index(1), &
               SPtranspose(Asub2intf(i,i)) ) ! 2 stands for right

       endif
    enddo

    !! filling the SUBDOMAIN-to-SUBDOMAIN connectivity
    do i=1, osm_opt%N_subd
       if (i > 1 ) then
          call insertblock2spmat(AaugOff, &
               osm_opt%subd_interf2Gdof(i,1)%index(1), &
               osm_opt%subd_interf2Gdof(i-1,2)%index(1), &
               (1.0d0 - gamma) * Aintf(i-1) )     ! 1 stands for left
          call insertblock2spmat(AaugOff, &
               osm_opt%subd_interf2Gdof(i,1)%index(1), &
               osm_opt%subd_primal2Gdof(i-1)%index(1), &
               SPtranspose(Asub2intf(i-1,i-1)) ) ! 1 stands for left
       endif

       if (i < osm_opt%N_subd) then
          !
          call insertblock2spmat(AaugOff, &
               osm_opt%subd_interf2Gdof(i,2)%index(1), &
               osm_opt%subd_interf2Gdof(i+1,1)%index(1), &
               (1.0d0 - gamma) * Aintf(i) )     ! 2 stands for right
          call insertblock2spmat(AaugOff, &
               osm_opt%subd_interf2Gdof(i,2)%index(1), &
               osm_opt%subd_primal2Gdof(i+1)%index(1), &
               SPtranspose(Asub2intf(i+1,i)) ) ! 2 stands for right
          !
       endif
       
    enddo

  call deleteblankspmat(AaugD)
  call deleteblankspmat(AaugOff)
  !! 
  end subroutine build_augmatrix

  
end module fd
