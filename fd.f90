!> fd is a module that contains subroutines for
!! building finite difference discretizations and
!! domain decomposition components.
module fd
  use geo2d
  use sparsemat

  implicit none

contains

  !> fd_stiffness builds stiffness matrix
  !! based on the triangulation Th and writes
  !! in Asp, a sparse matrix.
  !! Asp should NOT had been initialized.
  subroutine fd_stiffness(Th,Asp)
    type(mesh2d), intent(in)    :: Th
    type(spmat),  intent(inout) :: Asp

    integer :: n_dof            ! total number of dofs
    real    :: hx               ! hx is the mesh size in x direction
    real    :: hy               ! hy is the mesh size in y direction

    ! in order to save time we compute -1/hx^2, +2/hx^2, -1/hy^2, +2/hy^2
    real    :: mHX, p2HX, mHY, p2HY

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
          call insert2spmat(res(i), j, geo_dd%subdm(i)%index(j), 1.0)
       end do
       call deleteblankspmat(res(i))
    end do

    !! building restriction for the interfaces
    do i=1, N_intf
       call buildspmat(resGAM(i), geo_dd%interf(i)%n_index, Th%n_vt, &
            geo_dd%interf(i)%n_index )
       do j=1, geo_dd%interf(i)%n_index
          call insert2spmat(resGAM(i), j, geo_dd%interf(i)%index(j), 1.0)
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

    real    :: hx, hy
    real    :: mHX, p2HX, mHY, p2HY

    real :: Eta

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

    Eta = 1.0/ sqrt(hx**2 + hy**2)

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
end module fd
