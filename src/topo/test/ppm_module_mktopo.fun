test_suite ppm_module_mktopo


#ifdef __MPI
    INCLUDE "mpif.h"
#endif

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: pi = 3.1415926535897931_mk
REAL(MK),PARAMETER              :: skin = 0.0_MK
INTEGER,PARAMETER               :: ndim=2
INTEGER,PARAMETER               :: pdim=2
INTEGER                         :: decomp,assig,tolexp
REAL(MK)                        :: tol
REAL(MK)                        :: cutoff = 0.05_mk
INTEGER                         :: info,comm,rank,nproc
INTEGER                         :: topoid
INTEGER                         :: np = 100000
INTEGER                         :: mp
INTEGER                         :: newnp
REAL(MK),DIMENSION(:,:),POINTER :: xp
REAL(MK),DIMENSION(:  ),POINTER :: rcp
REAL(MK),DIMENSION(:,:),POINTER :: wp
REAL(MK),DIMENSION(:  ),POINTER :: min_phys,max_phys,h,p_h
REAL(MK),DIMENSION(:  ),POINTER :: len_phys
REAL(MK),DIMENSION(:  ),POINTER :: ghostlayer
INTEGER, DIMENSION(:  ),POINTER :: ghostsize
INTEGER                         :: i,j,k,sum1,sum2
INTEGER                         :: p_i
INTEGER, DIMENSION(6)           :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost
INTEGER, DIMENSION(:  ),POINTER :: nm
INTEGER                         :: seedsize
INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
REAL(MK), DIMENSION(:),ALLOCATABLE :: randnb
INTEGER                          :: isymm = 0
LOGICAL                          :: lsymm = .FALSE.,ok
REAL(MK)                         :: t0,t1,t2,t3

    init

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_init

        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),ghostlayer(2*ndim),&
            &         nm(ndim),h(ndim),p_h(ndim),STAT=info)

        min_phys(1:ndim) = 0.0_MK
        max_phys(1:ndim) = 1.0_MK
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = cutoff
        bcdef(1:6) = ppm_param_bcdef_periodic

        NULLIFY(xp,rcp,wp)

#ifdef __MPI
        comm = MPI_COMM_WORLD
        CALL MPI_Comm_rank(comm,rank,info)
        CALL MPI_Comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = INT(LOG10(EPSILON(1.0_MK)))+10
        CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        USE ppm_module_finalize

        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize


    setup

        CALL RANDOM_SEED(size=seedsize)
        ALLOCATE(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        CALL RANDOM_SEED(put=seed)
        ALLOCATE(randnb((1+ndim)*np),STAT=info)
        CALL RANDOM_NUMBER(randnb)
        ALLOCATE(xp(ndim,np),rcp(np),wp(pdim,np),STAT=info)

    end setup


    teardown

        DEALLOCATE(xp,rcp,wp,STAT=info)
        DEALLOCATE(seed,randnb)

    end teardown

    test cuboid
        ! test cuboid decomposition

        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_write
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo

        start_subroutine("cuboid")

        !----------------
        ! create particles
        !----------------
        CALL RANDOM_NUMBER(xp)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_bisection
        !decomp = ppm_param_decomp_xpencil
        !assig  = ppm_param_assign_internal
        assig  = ppm_param_assign_metis_comm

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)

        Assert_Equal(info,0)

        stdout('ppm_topo(topoid)%t%nsubs','ppm_topo(topoid)%t%nsublist')

        end_subroutine()

    end test

    test xp_null_cuboid
      ! ppm trac ticket #93 (https://ppm.inf.ethz.ch/trac/ticket/93)
      USE ppm_module_util_dbg

      INTEGER :: meshid = -1

      DEALLOCATE(xp)
      xp => null()

      nm = 32 * nproc

      CALL ppm_mktopo(topoid, meshid,     &
                      xp, 0,              &
                      decomp, assig,      &
                      min_phys, max_phys, &
                      bcdef, ghostsize,  &
                      cost, Nm, info)

      Assert_Equal(info, 0)

!       CALL ppm_dbg_print(topoid,0.0_MK,1,1,info)

    end test


end test_suite
