test_suite ppm_module_util_netstat



#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol
integer                         :: info,comm,rank,nproc
integer                         :: topoid
integer                         :: np = 1000
real(mk),dimension(:,:),pointer :: xp => NULL()
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(mk),dimension(:  ),pointer :: len_phys => NULL()
integer                         :: i,j,k
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
real(mk)                        :: latency,bandwidth
logical                         :: ok

    init

        use ppm_module_data
        use ppm_module_init
        use ppm_module_topo_typedef

        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         stat=info)

        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        tol = epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))

        nullify(xp)

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys)

    end finalize


    setup

        allocate(xp(ndim,np),stat=info)

    end setup


    teardown

        deallocate(xp,stat=info)

    end teardown

    test netstat
        ! test netstat

        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_topo_check
        use ppm_module_test

        !----------------
        ! create particles
        !----------------

        call part_init(xp,np,min_phys,max_phys,info)


        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               0.1_mk,cost,info)

        call ppm_netstat(topoid,latency,bandwidth,info)
        write(*,'(A,I3,A,F12.6,A,F10.2,A)') 'rank: ',ppm_rank,&
        &                         '  latency: ',latency*1000.0_mk,&
        &                         ' ms  bandwidth: ',bandwidth/REAL(1024**2,mk),&
        &                         ' MB/s'
        if (info.ge.-2) then
            ok = .true.
        endif
        assert_true(ok)

    end test


end test_suite
