test_suite ppm_module_util_netstat

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
INTEGER, PARAMETER               :: ndim=2
INTEGER                         :: decomp,assig,tolexp
REAL(MK)                        :: tol
INTEGER                         :: info,comm,rank,nproc
INTEGER                         :: topoid
INTEGER                         :: np = 1000
REAL(MK),DIMENSION(:,:),POINTER :: xp => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: len_phys => NULL()
INTEGER                         :: i,j,k
INTEGER, DIMENSION(6)           :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
REAL(MK)                        :: latency,bandwidth
LOGICAL                         :: ok

    init

        USE ppm_module_data
        USE ppm_module_init
        USE ppm_module_topo_typedef

        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),STAT=info)

        min_phys(1:ndim) = 0.0_MK
        max_phys(1:ndim) = 1.0_MK
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        tol = EPSILON(1.0_MK)
        tolexp = INT(LOG10(EPSILON(1.0_MK)))

        NULLIFY(xp)

#ifdef __MPI
        comm = MPI_COMM_WORLD
        CALL MPI_Comm_rank(comm,rank,info)
        CALL MPI_Comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        USE ppm_module_finalize

        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys,len_phys)

    end finalize


    setup

        ALLOCATE(xp(ndim,np),STAT=info)

    end setup


    teardown

        DEALLOCATE(xp,STAT=info)

    end teardown

    test netstat
        ! test netstat

        USE ppm_module_data
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_test

        !----------------
        ! create particles
        !----------------

        CALL part_init(xp,np,min_phys,max_phys,info)


        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef,0.1_MK,cost,info)

        CALL ppm_netstat(topoid,latency,bandwidth,info)
        write(*,'(A,I3,A,F12.6,A,F10.2,A)') 'rank: ',ppm_rank,'  latency: ',&
        & latency*1000.0_MK,' ms  bandwidth: ',bandwidth/REAL(1024**2,mk),' MB/s'
        IF (info.GE.-2) then
            ok = .TRUE.
        ENDIF
        Assert_True(ok)

    end test


end test_suite
