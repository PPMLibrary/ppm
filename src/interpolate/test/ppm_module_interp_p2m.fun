test_suite ppm_module_interp_p2m



#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    integer, parameter              :: debug = 0
    integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
#ifdef __MPI
    integer, parameter              :: comm = mpi_comm_world
#endif
    integer                         :: ndim,nspec
    integer                         :: rank
    integer                         :: nproc
    integer                         :: decomp
    integer                         :: assig
    integer                         :: tolexp
    real(mk)                        :: tol
    integer                         :: info
    integer                         :: topoid,meshid
    real(mk),dimension(:,:),pointer :: xp,wp
    real(mk),dimension(:  ),pointer :: min_phys,max_phys,h
    integer, dimension(:  ),pointer :: ghostsize
    integer                         :: i,j,k,p_i,ai,aj,it,isub
    integer, dimension(6)           :: bcdef
    real(mk),dimension(:  ),pointer :: cost
    integer, dimension(:,:),pointer :: istart,ndata
    integer, dimension(:  ),pointer :: nm
    integer                         :: np,mp
    integer                         :: kernel
    real(mk),dimension(:,:,:,:  ), pointer :: field_wp2 ! field_wp(ldn,i,j,isub)
    real(mk),dimension(:,:,:,:,:), pointer :: field_wp3 ! field_wp(ldn,i,j,k,isub)
    real(mk),dimension(:  ),pointer :: field_x
    real(mk)                        :: maxm3
    integer                         :: seedsize
    integer,  dimension(:),allocatable :: seed
    !---- The following variables are needed for testing the moments
    integer, parameter              :: nmom2 = 10
    integer, parameter              :: nmom3 = 19
    integer, dimension(2,nmom2)     :: alpha2
    integer, dimension(3,nmom3)     :: alpha3
    real(mk),dimension(nmom2)       :: f_moments2, f_mom_global2, p_moments2
    real(mk),dimension(nmom3)       :: f_moments3, f_mom_global3, p_moments3

    data ((alpha2(ai,aj), ai=1,2), aj=1,nmom2) /0,0, 1,0, 0,1, 2,0, 0,2, &
    &                                           1,1, 3,0, 0,3, 2,1, 1,2/
    data ((alpha3(ai,aj), ai=1,3), aj=1,nmom3)               /0,0,0, &
    &                                           1,0,0, 0,1,0, 0,0,1, &
    &                      2,0,0, 0,2,0, 0,0,2, 1,1,0, 0,1,1, 1,0,1, &
    & 3,0,0, 0,3,0, 0,0,3, 0,1,2, 0,2,1, 1,2,0, 1,0,2, 2,1,0, 2,0,1/
!-------------------------- init testsuit -------------------------------------
    init
        use ppm_module_data

        tol = 100.0_mk*epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))
        bcdef(1:6) = ppm_param_bcdef_freespace

        allocate(min_phys(3),max_phys(3),ghostsize(3),&
        &         nm(3),h(3),field_x(3),stat=info)


#ifdef __MPI
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        nproc = 1
        rank = 0
#endif

    end init
!------------------------------------------------------------------------------


!------------------------- finalize testsuit ----------------------------------
    finalize

        deallocate(min_phys,max_phys,ghostsize,nm)

    end finalize
!------------------------------------------------------------------------------


!------------------------------ test setup ------------------------------------
    setup
        
        nullify(xp)
        nullify(wp)
        nullify(field_wp2)
        nullify(field_wp3)
        nullify(cost)

        np = 400*nproc
        mp = 0
        
        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

    end setup
!------------------------------------------------------------------------------
        

!--------------------------- test teardown ------------------------------------
    teardown
        
        use ppm_module_finalize

        call ppm_finalize(info)
        
        deallocate(xp,wp,stat=info)
        deallocate(field_wp2,field_wp3,stat=info)
        deallocate(seed)
        deallocate(cost)

    end teardown
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    test p2m_mp4_2d
        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_topo_get
        use ppm_module_interp_p2m
        use ppm_module_init
        use ppm_module_finalize
        use ppm_module_map

        implicit none
        integer, dimension(2)           :: maxndata
        INTEGER, DIMENSION(:  ), POINTER:: isublist 
        integer                         :: nsublist
        ndim = 2
        nspec = 1
        kernel = ppm_param_rmsh_kernel_mp4
        do i=1,ndim
            min_phys(i) = 0.0_mk
            max_phys(i) = 1.0_mk
            ghostsize(i) = 2
        enddo
        
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)
    
        allocate(xp(ndim,np),wp(nspec,np),stat=info)
        call random_number(xp)
        wp = 0.0_mk


        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0
        meshid = -1



        allocate(nm(ndim),stat=info)
        do i=1,ndim
            nm(i) = 32*nproc
        enddo

        call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,    &
        &               bcdef,ghostsize,cost,nm,info)
        
        call ppm_topo_get_meshinfo(topoid,meshid,nm,istart,ndata,maxndata,&
                        isublist,nsublist,info)


        allocate(field_wp2(nspec,(1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)),nsublist),&
        &        stat=info) ! 2d

        do i=1,ndim
            h(i) = (max_phys(i) - min_phys(i)) / real(ndata(i,1)-1,mk)
        enddo

        call ppm_map_part_global(topoid,xp,np,info) ! positions
        call ppm_map_part_push(wp,nspec,np,info)    ! strengths
        call ppm_map_part_send(np,mp,info)          ! send
        call ppm_map_part_pop(wp,nspec,np,mp,info)  ! strengths
        call ppm_map_part_pop(xp,ndim,np,mp,info)   ! positions
        np = mp

!        maxm3 = 0.0_mk
        do p_i=1,np
            !----------------
            ! p --> m
            !----------------
            wp(1,p_i) = 1.0_mk
            call ppm_interp_p2m(topoid,meshid,xp,np,wp,1,kernel,ghostsize,&
            &                   field_wp2,info)

            !----------------
            ! test p --> m
            !----------------
            f_moments2 = 0.0_mk
            p_moments2 = 0.0_mk
            do it=1,nsublist
            isub = isublist(it)
            do j = 1-ghostsize(2), ndata(2,isub)+ghostsize(2)
                do i = 1-ghostsize(1),ndata(1,isub)+ghostsize(1)
                    field_x(1) = min_phys(1) + h(1)*real(i-1,mk)
                    field_x(2) = min_phys(2) + h(2)*real(j-1,mk)
                    do aj = 1,nmom2
                        f_moments2(aj) = f_moments2(aj) + field_wp2(1,i,j,isub)* &
                        &         field_x(1)**alpha2(1,aj)*field_x(2)**alpha2(2,aj)
                    enddo
                enddo
            enddo
            enddo
#ifdef __MPI
            call mpi_allreduce(f_moments2,f_mom_global2,nmom2,ppm_mpi_kind,MPI_SUM,&
            &                  comm,info)
#else
            f_mom_global2 = f_moments2
#endif
            do aj = 1,nmom2
                p_moments2(aj) = xp(1,p_i)**alpha2(1,aj)*xp(2,p_i)**alpha2(2,aj)
            enddo
            do aj = 1,6
                assert_equal_within(f_mom_global2(aj),p_moments2(aj),tol)
!                if (abs(f_moments(aj) - p_moments(aj)) .GT. tol) then
!                    print *, 'particle pos:',     xp(:,p_i)
!                    print *, 'failed at moment: ', aj
!                    print *, 'field moments: ',   f_moments
!                    print *, 'particle moments: ',p_moments
!                    stop 'ERROR: p2m interpolation: moments not conserved.'
!                endif
            enddo
!            do aj = 7,10
!                if (abs(f_moments(aj) - p_moments(aj)) .GT. maxm3) then
!                    maxm3 = abs(f_moments(aj) - p_moments(aj))
!                endif
!            enddo
            wp(1,p_i) = 0.0_mk
        enddo

!        print *, 'Maximum 3rd moment difference / h^3', maxm3/h**3
        
    end test
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    test p2m_mp4_3d
        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_topo_get
        use ppm_module_interp_p2m
        use ppm_module_init
        use ppm_module_finalize
        use ppm_module_map

        implicit none
        integer, dimension(3)           :: maxndata
        INTEGER, DIMENSION(:  ), POINTER:: isublist 
        integer                         :: nsublist
        
        ndim = 3
        nspec = 1
        kernel = ppm_param_rmsh_kernel_mp4
        do i=1,ndim
            min_phys(i) = 0.0_mk
            max_phys(i) = 1.0_mk
            ghostsize(i) = 2
        enddo
        
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)
    
        allocate(xp(ndim,np),wp(nspec,np),stat=info)
        call random_number(xp)
        wp = 0.0_mk


        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0
        meshid = -1



        allocate(nm(ndim),stat=info)
        do i=1,ndim
            nm(i) = 32*nproc
        enddo

        call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,    &
        &               bcdef,ghostsize,cost,nm,info)
        
        call ppm_topo_get_meshinfo(topoid,meshid,istart,ndata,maxndata,&
                        isublist,nsublist,info)
        

        allocate(field_wp3(nspec,(1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)), &
        &        (1-ghostsize(3)):(maxndata(3)+ghostsize(3)),nsublist),&
        &       stat=info) ! 3d

        do i=1,ndim
            h(i) = (max_phys(i) - min_phys(i)) / real(ndata(i,1)-1,mk)
        enddo

        call ppm_map_part_global(topoid,xp,np,info) ! positions
        call ppm_map_part_push(wp,nspec,np,info)    ! strengths
        call ppm_map_part_send(np,mp,info)          ! send
        call ppm_map_part_pop(wp,nspec,np,mp,info)  ! strengths
        call ppm_map_part_pop(xp,ndim,np,mp,info)   ! positions
        np = mp

!        maxm3 = 0.0_mk
        do p_i=1,np
            !----------------
            ! p --> m
            !----------------
            wp(1,p_i) = 1.0_mk
            call ppm_interp_p2m(topoid,meshid,xp,np,wp,1,kernel,ghostsize,&
            &                   field_wp3,info)

            !----------------
            ! test p --> m
            !----------------
            f_moments3 = 0.0_mk
            p_moments3 = 0.0_mk
            do it=1,nsublist
            isub = isublist(it)
            do k = 1-ghostsize(3), ndata(3,isub)+ghostsize(3)
                do j = 1-ghostsize(2), ndata(2,isub)+ghostsize(2)
                    do i = 1-ghostsize(1),ndata(1,isub)+ghostsize(1)
                        field_x(1) = min_phys(1) + h(1)*real(i-1,mk)
                        field_x(2) = min_phys(2) + h(2)*real(j-1,mk)
                        field_x(3) = min_phys(3) + h(3)*real(k-1,mk)
                        do aj = 1,nmom3
                            f_moments3(aj) = f_moments3(aj) + &
                                field_wp3(1,i,j,k,isub)* &
                            &   field_x(1)**alpha3(1,aj)*&
                            &   field_x(2)**alpha3(2,aj)*&
                            &   field_x(3)**alpha3(3,aj)
                        enddo
                    enddo
                enddo
            enddo
            enddo
#ifdef __MPI
            call mpi_allreduce(f_moments3,f_mom_global3,nmom3,ppm_mpi_kind,MPI_SUM,&
            &                  comm,info)
#else
            f_mom_global3 = f_moments3
#endif
            do aj = 1,nmom3
                p_moments3(aj) = xp(1,p_i)**alpha3(1,aj)*&
                &                xp(2,p_i)**alpha3(2,aj)*&
                &                xp(3,p_i)**alpha3(3,aj)
            enddo
            do aj = 1,10
                assert_equal_within(f_mom_global3(aj),p_moments3(aj),tol)
!                if (abs(f_moments(aj) - p_moments(aj)) .GT. tol) then
!                    print *, 'particle pos:',     xp(:,p_i)
!                    print *, 'failed at moment: ', aj
!                    print *, 'field moments: ',   f_moments
!                    print *, 'particle moments: ',p_moments
!                    stop 'ERROR: p2m interpolation: moments not conserved.'
!                endif
            enddo
!            do aj = 7,10
!                if (abs(f_moments(aj) - p_moments(aj)) .GT. maxm3) then
!                    maxm3 = abs(f_moments(aj) - p_moments(aj))
!                endif
!            enddo
            wp(1,p_i) = 0.0_mk
        enddo

!        print *, 'Maximum 3rd moment difference / h^3', maxm3/h**3
        
    end test
!------------------------------------------------------------------------------


end test_suite
