test_suite ppm_module_mktopo


#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer,parameter               :: pdim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol
real(mk)                        :: cutoff = 0.05_mk
integer                         :: info,comm,rank,nproc
integer                         :: topoid
integer                         :: np = 100000
integer                         :: mp
integer                         :: newnp
real(mk),dimension(:,:),pointer :: xp
real(mk),dimension(:  ),pointer :: rcp
real(mk),dimension(:,:),pointer :: wp
real(mk),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: ghostlayer
integer, dimension(:  ),pointer :: ghostsize
integer                         :: i,j,k,sum1,sum2
integer                         :: p_i
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
integer, dimension(:  ),pointer :: nm
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok
real(mk)                         :: t0,t1,t2,t3

    init

        use ppm_module_data
        use ppm_module_topo_typedef
        use ppm_module_init

        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),ghostlayer(2*ndim),&
            &         nm(ndim),h(ndim),p_h(ndim),stat=info)

        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = cutoff
        bcdef(1:6) = ppm_param_bcdef_periodic

        nullify(xp,rcp,wp)

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = INT(LOG10(EPSILON(1._mk)))+10
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize


    setup

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)
        allocate(randnb((1+ndim)*np),stat=info)
        call random_number(randnb)
        allocate(xp(ndim,np),rcp(np),wp(pdim,np),stat=info)

    end setup


    teardown

        deallocate(xp,rcp,wp,stat=info)
        deallocate(seed,randnb)

    end teardown

    test cuboid
        ! test cuboid decomposition

        use ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_write
        use ppm_module_topo_typedef
        use ppm_module_mktopo

        start_subroutine("cuboid")

        !----------------
        ! create particles
        !----------------
        call random_number(xp)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_bisection
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)

        Assert_Equal(info,0)

        stdout('ppm_topo(topoid)%t%nsubs','ppm_topo(topoid)%t%nsublist')

        end_subroutine()

    end test

    test xp_null_cuboid
      ! ppm trac ticket #93 (https://ppm.inf.ethz.ch/trac/ticket/93)
      use ppm_module_util_dbg

      integer :: meshid = -1

      deallocate(xp)
      xp => null()

      nm = 32 * nproc

      call ppm_mktopo(topoid, meshid,     &
                      xp, 0,              &
                      decomp, assig,      &
                      min_phys, max_phys, &
                      bcdef, ghostsize,  &
                      cost, Nm, info)

      Assert_Equal(info, 0)

!       call ppm_dbg_print(topoid,0.0_mk,1,1,info)

    end test


end test_suite
