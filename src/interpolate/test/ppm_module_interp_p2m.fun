test_suite ppm_module_interp_p2m


#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    ! copied from ppm_param.h
    INTEGER, PARAMETER :: ppm_param_device_cpu    = 1
    INTEGER, PARAMETER :: ppm_param_device_gpu    = 2
    INTEGER, PARAMETER :: ppm_param_rmsh_kernel_bsp2    = 1
    INTEGER, PARAMETER :: ppm_param_rmsh_kernel_mp4     = 2
    integer, parameter              :: debug = 0
    !integer, parameter              :: mk = kind(1.0d0) ! double-precision
    integer, parameter              ::  mk = kind(1.0e0) ! single-precision
#ifdef __MPI
    integer, parameter              :: comm = mpi_comm_world
#endif
    integer         , parameter     :: nspec=1
    real(mk)           ,parameter   :: pi=acos(-1.0_mk)
    integer                         :: rank
    integer                         :: nproc
    integer                         :: decomp
    integer                         :: assig
    integer                         :: tolexp
    real(mk)                        :: tol
    integer                         :: info
    integer                         :: ngp
    integer                         :: topoid,meshid
    real(mk),dimension(:,:),pointer :: xp => NULL()
    real(mk),dimension(:,:),pointer :: wp => NULL()
    real(mk),dimension(:  ),pointer :: h => NULL()
    real(mk),dimension(:  ),pointer :: min_phys => NULL()
    real(mk),dimension(:  ),pointer :: max_phys => NULL()
    integer, dimension(:  ),pointer :: ghostsize => NULL()
    integer                         :: i,j,k,p_i,ai,aj,it,isub
    integer, dimension(6)           :: bcdef
    real(mk),dimension(:  ),pointer :: cost => NULL()
    integer, dimension(:,:),pointer :: istart => NULL()
    integer, dimension(:,:),pointer :: ndata => NULL()
    integer, dimension(:  ),pointer :: nm => NULL()
    integer                         :: np,mp
    integer                         :: kernel
    integer                         :: compdev
    real(mk),dimension(:,:,:,:  ), pointer :: field_wp2 => NULL() ! field_wp(ldn,i,j,isub)
    real(mk),dimension(:,:,:,:,:), pointer :: field_wp3 => NULL() ! field_wp(ldn,i,j,k,isub)
    real(mk),dimension(:  ),pointer :: field_x => NULL()
    real(mk)                        :: maxm3
    real(mk)                        :: L2_norm,Linf_norm
    real(mk)                        :: f2_norm,finf_norm
    real(mk)                        :: ana_sol,diff
    real(mk)                        :: var,mu
    real(mk)                        :: x_coor,y_coor,z_coor
    integer  :: seedsize
    integer,dimension(:), pointer   :: seed
    real(mk)                        :: roff
    integer                         :: testid


    real(mk) :: t_start, t_end, t_total_gpu, t_total_cpu
    integer  :: ntrial = 1

    real(mk) :: rn 
    integer  :: p,id
!-------------------------- init testsuit -------------------------------------
    init
        use ppm_module_data

        tol = 100.0_mk*epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))
        bcdef(1:6) = ppm_param_bcdef_freespace
        testid = 0
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

        np = 100000*nproc
        mp = 0
        
        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+42)
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

    end teardown
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    test p2m_2d({ngp: [17,33,65,129,257,513,1025], kernel: [ppm_param_rmsh_kernel_bsp2,ppm_param_rmsh_kernel_mp4], compdev: [ppm_param_device_cpu,ppm_param_device_gpu]})
        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_topo_get
        use ppm_module_interp
        use ppm_module_init
        use ppm_module_finalize
        use ppm_module_map
        use ppm_module_interp_p2m

        implicit none
        integer         , parameter     :: ndim=2
        integer, dimension(ndim)        :: midx
        integer, dimension(ndim)        :: maxndata
        integer, dimension(:  ), pointer:: isublist => NULL()
        integer                         :: nsublist
        integer                         :: navg
        real(mk)                        :: t,hmax
        
        do i=1,ndim
            min_phys(i) = 0.0_mk
            max_phys(i) = 1.0_mk
            ghostsize(i) = 2
        enddo
        mu  = (max_phys(1) - min_phys(1))/2.0_mk + min_phys(1)
        var = (max_phys(1) - min_phys(1))/15.0_mk
       

        call ppm_init(ndim,mk,tolexp,0,debug,info,99)
    
        allocate(nm(ndim),stat=info)
        do i=1,ndim
            nm(i) = ngp*nproc
        enddo

        np = (nm(1)-1)*(nm(2)-1)
        allocate(xp(ndim,np),wp(nspec,np),stat=info)

        do i=1,ndim
            h(i) = (max_phys(i) - min_phys(i)) / real(nm(i)-1,mk)
        enddo

        hmax = (max_phys(1) - min_phys(1)) / real(1024,mk)
        t = 0.1_mk*( log10(1.0_mk+hmax) )

        p = 0
        DO j = 1, nm(2) - 1
            midx(2) = j
            DO i = 1, nm(1) - 1
                midx(1) = i
                p = p + 1
                !DO id=1,ndim
                    !xp(id,p) = (REAL(midx(id),mk)-0.3_mk)*h(id) + &
                    !&          + min_phys(id)
                !ENDDO
                xp(1,p) = ((REAL(i,mk)-0.4_mk)*h(1) + min_phys(1)) * exp(t)
                xp(2,p) = ((REAL(j,mk)-0.4_mk)*h(2) + min_phys(2)) * exp(-t)
                wp(1,p) = EXP(-((xp(1,p)-mu)**2)/var)*EXP(-((xp(2,p)-mu)**2)/var)
            ENDDO
        ENDDO
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0
        meshid = -1
        testid = testid + 1
        print *, '-------------------------------------'
        print *,'test:',testid
        if (compdev.eq.ppm_param_device_cpu) then
            if (kernel.eq.ppm_param_rmsh_kernel_bsp2) then
                print *, 'P2M, 2D, BSP2, CPU', nm
            else if (kernel.eq.ppm_param_rmsh_kernel_mp4) then
                print *, 'P2M, 2D, MP4, CPU', nm
            end if
        else if (compdev.eq.ppm_param_device_gpu) then
            if (kernel.eq.ppm_param_rmsh_kernel_bsp2) then
                print *, 'P2M, 2D, BSP2, GPU', nm
            else if (kernel.eq.ppm_param_rmsh_kernel_mp4) then
                print *, 'P2M, 2D, MP4, GPU', nm
            endif
        endif

        call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,    &
        &               bcdef,ghostsize,cost,nm,info)
        
        call ppm_topo_get_meshinfo(topoid,meshid,nm,istart,ndata,maxndata,&
        &               isublist,nsublist,info)


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

        call ppm_interp_p2m(topoid,meshid,xp,np,wp,1,kernel,ghostsize,&
        &                   field_wp2,info, device=compdev)

        p = 0
        L2_norm   = 0.0_mk
        Linf_norm = 0.0_mk
        f2_norm = 0.0_mk
        finf_norm = 0.0_mk
        navg = 0
        DO j = 2, nm(2)-3
            DO i = 2, nm(1)-3
                navg = navg + 1
                x_coor  = REAL(i,mk)*h(1) + min_phys(1)
                y_coor  = REAL(j,mk)*h(2) + min_phys(2)
                ana_sol = EXP(-((x_coor-mu)**2)/var)*EXP(-((y_coor-mu)**2)/var)
                diff    = ABS(ana_sol - field_wp2(1,i+1,j+1,1))
                L2_norm = L2_norm + diff**2
                write(100+testid,'(3E14.6)'),x_coor,y_coor,diff
                IF(diff .GT. Linf_norm) THEN
                    Linf_norm = diff
                ENDIF
                f2_norm = f2_norm + ana_sol**2
                IF(ana_sol .GT. finf_norm) THEN
                    finf_norm = ana_sol
                ENDIF
            ENDDO
        ENDDO

        L2_norm = SQRT(L2_norm/navg)
        f2_norm = SQRT(f2_norm/navg)
        write(*,'(A,3E14.6)')'h,L2,Linf',h(1),L2_norm/f2_norm,Linf_norm/finf_norm


    end test
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    test p2m_3d({ngp: [9,17,33,65,129,257], kernel: [ppm_param_rmsh_kernel_bsp2,ppm_param_rmsh_kernel_mp4], compdev: [ppm_param_device_cpu,ppm_param_device_gpu]})
        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_topo_get
        use ppm_module_interp
        use ppm_module_init
        use ppm_module_finalize
        use ppm_module_map
        use ppm_module_interp_p2m

        implicit none
        integer         , parameter     :: ndim=3
        integer, dimension(ndim)        :: midx
        integer, dimension(ndim)        :: maxndata
        integer, dimension(:  ), pointer:: isublist => NULL()
        integer                         :: nsublist
        integer                         :: navg
        real(mk)                        :: t,hmax
        
        do i=1,ndim
            min_phys(i) = 0.0_mk
            max_phys(i) = 1.0_mk
            ghostsize(i) = 2
        enddo
        mu  = (max_phys(1) - min_phys(1))/2.0_mk + min_phys(1)
        var = (max_phys(1) - min_phys(1))/15.0_mk
        

        call ppm_init(ndim,mk,tolexp,0,debug,info,99)
    
        allocate(nm(ndim),stat=info)
        do i=1,ndim
            nm(i) = ngp*nproc
        enddo

        np = (nm(1)-1)*(nm(2)-1)*(nm(3)-1)
        allocate(xp(ndim,np),wp(nspec,np),stat=info)

        do i=1,ndim
            h(i) = (max_phys(i) - min_phys(i)) / real(nm(i)-1,mk)
        enddo
        
        hmax = (max_phys(1) - min_phys(1)) / real(256,mk)
        t = 0.25*(log10(1.0_mk+hmax))

        p = 0
        DO k = 1, nm(3) - 1
            midx(3) = k
            DO j = 1, nm(2) - 1
                midx(2) = j
                DO i = 1, nm(1) - 1
                    midx(1) = i
                    p = p + 1
                    !DO id=1,ndim
                        !xp(id,p) = (REAL(midx(id),mk)-0.3_mk)*h(id) + &
                        !&          + min_phys(id)
                    !ENDDO
                    xp(1,p) = ((REAL(i,mk)-0.5_mk)*h(1) + min_phys(1)) * exp(t)
                    xp(2,p) = ((REAL(j,mk)-0.5_mk)*h(2) + min_phys(2)) * exp(-t)
                    xp(3,p) = (REAL(k,mk)-0.3_mk)*h(3) + min_phys(3) 
                    wp(1,p) = EXP(-((xp(1,p)-mu)**2)/var)*&
                    &         EXP(-((xp(2,p)-mu)**2)/var)*&
                    &         EXP(-((xp(3,p)-mu)**2)/var)
                ENDDO
            ENDDO
        ENDDO
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0
        meshid = -1
        
        print *, '-------------------------------------'
        if (compdev.eq.ppm_param_device_cpu) then
            if (kernel.eq.ppm_param_rmsh_kernel_bsp2) then
                print *, 'P2M, 3D, BSP2, CPU', nm
            else if (kernel.eq.ppm_param_rmsh_kernel_mp4) then
                print *, 'P2M, 3D, MP4, CPU', nm
            end if
        else if (compdev.eq.ppm_param_device_gpu) then
            if (kernel.eq.ppm_param_rmsh_kernel_bsp2) then
                print *, 'P2M, 3D, BSP2, GPU', nm
            else if (kernel.eq.ppm_param_rmsh_kernel_mp4) then
                print *, 'P2M, 3D, MP4, GPU', nm
            endif
        endif

        call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,    &
        &               bcdef,ghostsize,cost,nm,info)
        
        call ppm_topo_get_meshinfo(topoid,meshid,nm,istart,ndata,maxndata,&
        &               isublist,nsublist,info)


        allocate(field_wp3(nspec,(1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)),                  &
        &        (1-ghostsize(3)):(maxndata(3)+ghostsize(3)),nsublist),        &
        &        stat=info) ! 3d

        do i=1,ndim
            h(i) = (max_phys(i) - min_phys(i)) / real(ndata(i,1)-1,mk)
        enddo

        call ppm_map_part_global(topoid,xp,np,info) ! positions
        call ppm_map_part_push(wp,nspec,np,info)    ! strengths
        call ppm_map_part_send(np,mp,info)          ! send
        call ppm_map_part_pop(wp,nspec,np,mp,info)  ! strengths
        call ppm_map_part_pop(xp,ndim,np,mp,info)   ! positions
        np = mp

        call ppm_interp_p2m(topoid,meshid,xp,np,wp,1,kernel,ghostsize,&
        &                   field_wp3,info, device=compdev)

        p = 0
        L2_norm   = 0.0_mk
        Linf_norm = 0.0_mk
        navg = 0
        DO k = 2, nm(3)-3
            DO j = 2, nm(2)-3
                DO i = 2, nm(1)-3
                    navg = navg + 1
                    x_coor  = REAL(i,mk)*h(1) + min_phys(1)
                    y_coor  = REAL(j,mk)*h(2) + min_phys(2)
                    z_coor  = REAL(k,mk)*h(3) + min_phys(3)
                    ana_sol = EXP(-((x_coor-mu)**2)/var)*&
                    &         EXP(-((y_coor-mu)**2)/var)*&
                    &         EXP(-((z_coor-mu)**2)/var)
                    diff    = ABS(ana_sol - field_wp3(1,i+1,j+1,k+1,1))
                    L2_norm = L2_norm + diff**2
                    IF(diff .GT. Linf_norm) THEN
                        Linf_norm = diff
                    ENDIF
                    f2_norm = f2_norm + ana_sol**2
                    IF(ana_sol .GT. finf_norm) THEN
                        finf_norm = ana_sol
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        L2_norm = SQRT(L2_norm/navg)
        f2_norm = SQRT(f2_norm/navg)
        write(*,'(A,3E14.6)')'h,L2,Linf',h(1),L2_norm/f2_norm,Linf_norm/finf_norm


    end test
!------------------------------------------------------------------------------



end test_suite
