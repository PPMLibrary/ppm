test_suite map_part

    init

        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),ghostlayer(2*ndim),&
            &         nm(ndim),h(ndim),p_h(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = max_rcp
        bcdef(1:6) = ppm_param_bcdef_periodic
        
        nullify(xp,rcp,wp)

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_init(info)
        call mpi_comm_rank(comm,rank,info)
#else
        rank = 0
#endif
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize

        call ppm_finalize(info)
        #ifdef __MPI
        call MPI_finalize(info)
        #endif

        deallocate(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize


    setup

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        allocate(randnb((1+ndim)*np),stat=info)
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)
        call random_number(randnb)

    end setup

    test global
        !----------------
        ! create particles
        !----------------

        allocate(xp(ndim,np),rcp(np),wp(pdim,np),stat=info)
        xp = 0.0_mk
        rcp = 0.0_mk

        !p_h = len_phys / real(npgrid,mk)
        !do j=1,npgrid
        !    do i=1,npgrid
        !        p_i = i + (j-1)*npgrid
        !        xp(1,p_i) = min_phys(1)+real(i-1,mk)*p_h(1)
        !        xp(2,p_i) = min_phys(2)+real(j-1,mk)*p_h(2)
        !        rcp(p_i) = min_rcp + (max_rcp-min_rcp)*randnb(p_i)
        !        do k=1,pdim
        !            wp(k,i) = rcp(i)*REAL(k,MK)
        !        enddo
        !    enddo
        !enddo
        do i=1,np
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)

        print *,rank,'global map:'
        call ppm_map_part_global(topoid,xp,np,info)
        call ppm_map_part_push(rcp,np,info)
        call ppm_map_part_push(wp,pdim,np,info)
        call ppm_map_part_send(np,newnp,info)
        call ppm_map_part_pop(wp,pdim,np,newnp,info)
        call ppm_map_part_pop(rcp,np,newnp,info)
        call ppm_map_part_pop(xp,ndim,np,newnp,info)
        np=newnp
        print *,rank,'done'

        call ppm_topo_check(topoid,xp,np,ok,info)

        assert_true(ok)

    end test




end test_suite
