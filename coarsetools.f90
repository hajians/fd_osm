!> coarsetools contains subroutines for handling coarse/fine
!! interactions.
module coarsetools
  use geo2d
  use sparsemat
  
  implicit none

  type coarseT
     !> initialization flag
     logical :: init = .false.

     !> fine/coarse factor in x- and y-directions
     integer :: FCfactor_x, FCfactor_y

     !> Coarse vertex to fine vertex
     integer, allocatable :: Cv2Fv(:)

     !> prolongation stencil
     double precision, allocatable :: prolongStencil(:,:)

     !> prolongation operator
     type(spmat) :: prolongation
     
     !> restriction operator
     type(spmat) :: restriction
     
  end type coarseT

  private coarseVtx2fineVtx
  private buildprolongStencil
  private buildprolongation
  
contains

  !> if fine and coarse meshes constructed in a nested fashion then
  !! coarseVtx2fineVtx provides vertices of coarse mesh that sit on
  !! the fine ones.
  subroutine coarseVtx2fineVtx(ThC,ThF,coarse)
    type(mesh2d),  intent(in)    :: ThC
    type(mesh2d),  intent(in)    :: ThF
    type(coarseT), intent(inout) :: coarse

    integer :: Cvtx, Fvtx
    integer :: iC, jC, idxF, jdxF
    
    allocate( coarse%Cv2Fv(ThC%n_vt) ) ! allocate the array

    coarse%FCfactor_x = ThF%nx / ThC%nx 
    coarse%FCfactor_y = ThF%ny / ThC%ny

    do jC=1, ThC%n_vty
       do iC=1, ThC%n_vtx

          Cvtx = iC + (jC-1)*ThC%n_vtx ! Coarse vertex index

          idxF = coarse%FCfactor_x * iC
          jdxF = coarse%FCfactor_y * jC

          Fvtx = idxF + (jdxF-1) * ThF%n_vtx

          coarse%Cv2Fv(Cvtx) = Fvtx    ! done
       enddo
    enddo
    
  end subroutine coarseVtx2fineVtx

  !> builds prolongation stencil
  subroutine buildprolongStencil(mycase,coarse)
    integer,       intent(in)    :: mycase
    type(coarseT), intent(inout) :: coarse

    select case (mycase)
    case (1)                    ! 9-points interpolation
       allocate( coarse%prolongStencil(3,3) )

       coarse%prolongStencil(1,1) = 0.25d0;       coarse%prolongStencil(3,3) = 0.25d0;
       coarse%prolongStencil(1,3) = 0.25d0;       coarse%prolongStencil(3,1) = 0.25d0;

       coarse%prolongStencil(1,2) = 0.50d0;       coarse%prolongStencil(2,1) = 0.50d0;
       coarse%prolongStencil(2,3) = 0.50d0;       coarse%prolongStencil(3,2) = 0.50d0;

       coarse%prolongStencil(2,2) = 1.00d0;
    case (2)
       allocate( coarse%prolongStencil(3,3) )

       coarse%prolongStencil(1,1) = 0.50d0;       coarse%prolongStencil(3,3) = 0.50d0;
       coarse%prolongStencil(1,3) = 0.00d0;       coarse%prolongStencil(3,1) = 0.00d0;

       coarse%prolongStencil(1,2) = 0.50d0;       coarse%prolongStencil(2,1) = 0.50d0;
       coarse%prolongStencil(2,3) = 0.50d0;       coarse%prolongStencil(3,2) = 0.50d0;

       coarse%prolongStencil(2,2) = 1.00d0;
       
    end select
    
  end subroutine buildprolongStencil

  !> builds prolongation operator based prolongation stencil
  subroutine buildprolongation(ThC, ThF, coarse)
    type(mesh2d),  intent(in)    :: ThC
    type(mesh2d),  intent(in)    :: ThF
    type(coarseT), intent(inout) :: coarse

    integer :: N_array
    integer :: Cvt, Fvt         ! Coarse vertex, fine vertex
    integer :: adjV, diagV      ! adjacent and diagonal vertex

    if ( allocated(coarse%prolongStencil).eqv..false. ) then
       write(6,*) "error in buildprolongation"
       stop
    endif

    N_array = size( coarse%prolongStencil ) * ThC%n_vt

    call buildspmat(coarse%prolongation, ThF%n_vt, ThC%n_vt, N_array)

    do Cvt=1, ThC%n_vt
       Fvt = coarse%Cv2Fv(Cvt)  ! finding corresponding fine vertex

       ! inserting the center of stencil
       call insert2spmat(coarse%prolongation, Fvt, Cvt, coarse%prolongStencil(2,2) )
       !! inserting adjacent vertices
       !! left
       adjV = ThF%v2v(Fvt,1)
       call insert2spmat(coarse%prolongation, adjV, Cvt, coarse%prolongStencil(2,1) )
       !! right
       adjV = ThF%v2v(Fvt,2)
       call insert2spmat(coarse%prolongation, adjV, Cvt, coarse%prolongStencil(2,3) )
       !! south
       adjV = ThF%v2v(Fvt,3)
       call insert2spmat(coarse%prolongation, adjV, Cvt, coarse%prolongStencil(1,2) )
       !! north
       adjV = ThF%v2v(Fvt,4)
       call insert2spmat(coarse%prolongation, adjV, Cvt, coarse%prolongStencil(3,2) )

       !! north-east
       diagV = ThF%v2v(Fvt,4) - 1
       call insert2spmat(coarse%prolongation, diagV, Cvt, coarse%prolongStencil(3,1) )
       !! north-west
       diagV = ThF%v2v(Fvt,4) + 1
       call insert2spmat(coarse%prolongation, diagV, Cvt, coarse%prolongStencil(3,3) )
       !! south-east
       diagV = ThF%v2v(Fvt,3) - 1
       call insert2spmat(coarse%prolongation, diagV, Cvt, coarse%prolongStencil(1,1) )
       !! south-west
       diagV = ThF%v2v(Fvt,3) + 1
       call insert2spmat(coarse%prolongation, diagV, Cvt, coarse%prolongStencil(1,3) )

    enddo

    call deleteblankspmat(coarse%prolongation)
    
  end subroutine buildprolongation
  
  !> builds coarse tools
  subroutine buildcoarsetools(ThC,ThF,coarse)
    type(mesh2d),  intent(in)    :: ThC
    type(mesh2d),  intent(in)    :: ThF
    type(coarseT), intent(inout) :: coarse

    integer :: prolong_id

    if ( (ThC%init.eqv..false.) .or. (ThF%init.eqv..false.) ) then
       write(6,*) "error in buildcoarsetools", ThC%init, ThF%init
       stop
    endif

    if ( coarse%init.eqv..true. ) then
       write(6,*) "error in buildcoarse: coarse is initialized"
       stop
    endif

    prolong_id = 2
    
    call coarseVtx2fineVtx(ThC,ThF,coarse)

    call buildprolongStencil(prolong_id,coarse)
    call buildprolongation(ThC,ThF,coarse)
    
  end subroutine buildcoarsetools
  
  
end module coarsetools
