
    ! !! building local matrices Asub(:)
    ! do i=1, N_subd
    !    sub_dof = geo_dd%subdm(i)%n_index
    !    call buildspmat(Asub(i), sub_dof, sub_dof, 5*sub_dof)
    !    ! filling the matrix Asub(i)
    !    do n=1, sub_dof
    !       do m=1, sub_dof
    !          call readspmat(Asp, geo_dd%subdm(i)%index(n), geo_dd%subdm(i)%index(m), &
    !               val, flag)
    !          if ( flag > 0 ) then
    !             call insert2spmat(Asub(i),n,m,val)
    !          endif
    !       enddo
    !    enddo
    !    call deleteblankspmat(Asub(i))
    !    ! end filling the matrix Asub(i)
    ! enddo

    ! !! building local matrices Asubd2intf(:)
    ! do intf=1, N_intf
    !    intf_dof = geo_dd%interf(intf)%n_index
    !    !
    !    do j=1, 2             ! two subdomains per interface
    !       subd = geo_dd%intf2subd(intf,j)
    !       sub_dof = geo_dd%subdm(subd)%n_index
    !       call buildspmat(Asub2intf(subd,intf), sub_dof, intf_dof, 5*sub_dof)
    !       ! filling the matrix Asub(i)
    !       !write(6,*) intf, subd
    !       do n=1, sub_dof
    !          do m=1, intf_dof
    !             ! 
    !             call readspmat(Asp, &
    !                  geo_dd%subdm(subd)%index(n), geo_dd%interf(intf)%index(m), &
    !                  val, flag)
    !             if ( flag > 0 ) then
    !                call insert2spmat(Asub2intf(subd,intf),n,m,val)
    !             endif
    !          enddo
    !       enddo
    !       call deleteblankspmat( Asub2intf(subd,intf) )
    !       ! end filling the matrix Asub(i)       
    !    end do
    ! end do
