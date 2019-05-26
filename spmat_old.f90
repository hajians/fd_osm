!> @file spmat.f90 contains subroutines that manipulate a sparse matrix
!> @author Soheil Hajian
!! last update: July 21, 2014

module sparsemat

  implicit none

  !> a parameter to determine how many digits we need
  integer, parameter :: spkind = selected_int_kind (8)

  !> spmat is a type for sparse matrices in i,j,k format.
  type spmat

     logical :: init = .false.

     !> size of the sparse matrix
     integer (kind=spkind) :: N_i, N_j

     !> size of the i, j, k arrays
     integer (kind=spkind) :: N_array

     !> number of non-zero elements of the i, j, k arrays
     integer (kind=spkind) :: N_nonzero

     integer (kind=spkind), allocatable :: i(:)
     integer (kind=spkind), allocatable :: j(:)

     real, allocatable :: k(:)

  end type spmat


  !! private subroutine
  
contains

  !> buildspmat initiates a sparse matrix of size N_i, N_j
  !! with possible size of the i, j, k arrays.
  subroutine buildspmat(Asp, N_i, N_j, N_array)
    type(spmat),       intent(inout) :: Asp
    integer (kind=spkind),           intent(in)    :: N_i
    integer (kind=spkind),           intent(in)    :: N_j
    integer (kind=spkind), optional, intent(in)    :: N_array

    integer (kind=spkind) :: dummy

    if ( Asp%init.eqv..true. ) then
       write(6,*) "error in buildspmat: Asp is init."
       stop
    endif

    Asp%N_i = N_i
    Asp%N_j = N_j

    Asp%N_nonzero = 0

    if ( present(N_array) ) then
       dummy = N_array
    else
       dummy = 1
    endif

    allocate( Asp%i( dummy ) )
    allocate( Asp%j( dummy ) )
    allocate( Asp%k( dummy ) )

    Asp%k = 0.0

    Asp%N_array = dummy

    Asp%init = .true.

  end subroutine buildspmat

  !> clearspmat clears a sparse matrix type
  subroutine clearspmat(Asp)
    type(spmat), intent(inout) :: Asp

    if (Asp%init.eqv..false.) then
       write(6,*) "error in clearspmat"
       stop
    endif

    deallocate( Asp%i, Asp%j, Asp%k )

    Asp%N_i       = 0
    Asp%N_j       = 0
    Asp%N_array   = 0
    Asp%N_nonzero = 0

    Asp%init = .false.
  end subroutine clearspmat

  !> write sp mat in COO format
  subroutine writespmat(filename,Asp)
    character(len=*), intent(in)    :: filename
    type(spmat),      intent(in)    :: Asp

    integer (kind=spkind) :: n
    integer, parameter :: MyUnit = 11

    if ( Asp%init.eqv..false. ) then
       write(6,*) "Sparse matrix is not init in showspmat"
       stop
    endif

    open(MyUnit, file=filename)
    !! write dimension of the matrix and # of non-zero elements
    write(MyUnit,*) Asp%N_i, Asp%N_j, Asp%N_nonzero

    do n=1, Asp%N_nonzero
       write(MyUnit,*) Asp%i(n), Asp%j(n), Asp%k(n)
    enddo

    close(MyUnit)

  end subroutine writespmat

  !> show sparse matrix
  subroutine showspmat(Asp)
    type(spmat), intent(in) :: Asp

    integer (kind=spkind) :: n

    if ( Asp%init.eqv..false. ) then
       write(6,*) "Sparse matrix is not init in showspmat"
       stop
    endif

    write(6,*) "size of sparse matrix:     ", Asp%N_i, "*", Asp%N_j
    write(6,*) "number of non-zero entries:", Asp%N_nonzero
    write(6,*) "size of i, j, k arrays:     ", Asp%N_array

    do n=1, Asp%N_nonzero
       write(6,*) Asp%i(n), ",", Asp%j(n), ",", Asp%k(n)
    enddo

  end subroutine showspmat

  !> delete blank space of the arrays
  subroutine deleteblankspmat(Asp)
    type(spmat), intent(inout) :: Asp

    integer (kind=spkind), dimension( Asp%N_nonzero ) :: dummy_int
    real,                  dimension( Asp%N_nonzero ) :: dummy_real

    dummy_int = Asp%i(1:Asp%N_nonzero)
    deallocate( Asp%i )
    allocate( Asp%i( Asp%N_nonzero ) )
    Asp%i = dummy_int

    dummy_int = Asp%j(1:Asp%N_nonzero)
    deallocate( Asp%j )
    allocate( Asp%j( Asp%N_nonzero ) )
    Asp%j = dummy_int

    dummy_real = Asp%k(1:Asp%N_nonzero)
    deallocate( Asp%k )
    allocate( Asp%k( Asp%N_nonzero ) )
    Asp%k = dummy_real  

    Asp%N_array = Asp%N_nonzero
  end subroutine deleteblankspmat

  !> insert blank space into the arrays
  subroutine insertblank2spmat(Asp,N_blank)
    type(spmat), intent(inout) :: Asp
    integer (kind=spkind) , intent(in) :: N_blank

    integer (kind=spkind), dimension( Asp%N_nonzero ) :: dummy_int
    real,                  dimension( Asp%N_nonzero ) :: dummy_real

    Asp%N_array = Asp%N_array + N_blank

    dummy_int = Asp%i(1:Asp%N_nonzero)
    !
    deallocate( Asp%i )
    allocate( Asp%i(Asp%N_array) )
    Asp%i(1:Asp%N_nonzero) = dummy_int

    dummy_int = Asp%j(1:Asp%N_nonzero)
    !
    deallocate( Asp%j )
    allocate( Asp%j(Asp%N_array) )
    Asp%j(1:Asp%N_nonzero) = dummy_int

    dummy_real = Asp%k(1:Asp%N_nonzero)
    !
    deallocate( Asp%k )
    allocate( Asp%k(Asp%N_array) )
    Asp%k(1:Asp%N_nonzero) = dummy_real

  end subroutine insertblank2spmat

  !> insert an element in the sparse matrix
  subroutine insert2spmat(Asp,i,j,val)
    type(spmat),           intent(inout) :: Asp
    integer (kind=spkind), intent(in)    :: i
    integer (kind=spkind), intent(in)    :: j
    real,                  intent(in)    :: val

    integer (kind=spkind) :: N_blank

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in insert2spmat: init"
       stop
    endif

    if ( (Asp%N_i < i) .or. (Asp%N_j < j) ) then
       write(6,*) "error in insert2spmat: N_i < i or N_j < j"
       stop
    endif

    if ( Asp%N_array.eq.Asp%N_nonzero ) then
       N_blank = max( Asp%N_nonzero/2, 1 )
       call insertblank2spmat(Asp,N_blank)
    endif

    Asp%N_nonzero = Asp%N_nonzero + 1

    Asp%i( Asp%N_nonzero ) = i
    Asp%j( Asp%N_nonzero ) = j
    Asp%k( Asp%N_nonzero ) = val

  end subroutine insert2spmat

  !> read an entry of the sparse matrix.
  !! we just do a linear search in the arrays.
  subroutine readspmat(Asp,i,j,val,flag)
    type(spmat),           intent(in)    :: Asp
    integer (kind=spkind), intent(in)    :: i
    integer (kind=spkind), intent(in)    :: j
    real,                  intent(inout) :: val
    integer,               intent(inout) :: flag

    integer :: n

    if ( (i>Asp%N_i).or.(j>Asp%N_j).or.(Asp%init.eqv..false.) ) then
       write(6,*) "error in readspmat: "
    endif

    val  = 0.0
    flag = 0
    do n=1, Asp%N_nonzero
       !
       if ( Asp%i(n).eq.i ) then
          if ( Asp%j(n).eq.j ) then
             !
             val = val + Asp%k(n)
             flag = flag + 1
             !
          endif
       endif
       !
    enddo

    if (flag.ge.2) then
       write(6,*) "warning in readspmat: multiple entries found."
    endif
  end subroutine readspmat

  !> copysp2sp copies a sparse matrix into another one.
  subroutine copysp2sp(Asp,Bsp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(inout) :: Bsp
    
    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in copysp2sp: ", Asp%init
       stop
    endif
    !
    if (Bsp%init.eqv..true. ) then
       write(6,*) "warning in copysp2sp: Bsp is erased and re-init."
       call clearspmat( Bsp )
    endif

    call buildspmat(Bsp, Asp%N_i, Asp%N_j, Asp%N_array)

    Bsp%N_nonzero = Asp%N_nonzero

    Bsp%i(1:Asp%N_nonzero) = Asp%i(1:Asp%N_nonzero)
    Bsp%j(1:Asp%N_nonzero) = Asp%j(1:Asp%N_nonzero)
    Bsp%k(1:Asp%N_nonzero) = Asp%k(1:Asp%N_nonzero)

  end subroutine copysp2sp

  !> sub_redundant finds multiple entries at the same position
  !! in the sparse matrix form and sum them. This is feature is used
  !! in finite element method, where local assembly is used.
  !! you can use same matrix for both arguments:
  !! call sum_redundant(A1,A1)
  subroutine sum_redundant(Asp,Bsp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(inout) :: Bsp

    integer,     allocatable, dimension(:)   :: flag
    type(spmat)                              :: Acopy

    integer :: m, n
    integer :: i, j
    real    :: val

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in sum_redundant: ", Asp%init
       stop
    endif

    call copysp2sp(Asp,Acopy)   ! copy Asp to Acopy and use Acopy from now on

    if (Bsp%init.eqv..true. ) then
       write(6,*) "warning in sum_redundant: Bsp is erased and re-init."
       call clearspmat( Bsp )
    endif

    call buildspmat(Bsp, Acopy%N_i, Acopy%N_j, Acopy%N_nonzero)

    allocate( flag(Acopy%N_nonzero) ) ! allocate flag
    flag = -1                         ! -1 means not visited

    do n=1, Acopy%N_nonzero
       !
       if ( flag(n).eq.-1 ) then ! if not visited
          i   = Acopy%i(n)
          j   = Acopy%j(n)
          val = Acopy%k(n)
          !
          do m=(n+1), Acopy%N_nonzero
             if ( (Acopy%i(m).eq.i) .and. (Acopy%j(m).eq.j) ) then
                flag(m) = 0
                val = val + Acopy%k(m)
             endif
          enddo
          !
          call insert2spmat(Bsp,i,j,val) ! finally insert into Bsp
       endif
       !
    enddo

  end subroutine sum_redundant
end module sparsemat
