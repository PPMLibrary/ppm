test_suite ppm_module_util_color

    USE ppm_module_data
    USE ppm_module_substart
    USE ppm_module_substop

    INTEGER, PARAMETER                    :: debug = 0
    INTEGER, PARAMETER                    :: MK = KIND(1.0d0)
    INTEGER, PARAMETER                    :: ndim=3

    REAL(MK)                              :: cutoff = 0.05_MK
    REAL(MK), DIMENSION(:,:), POINTER     :: xp
    REAL(MK), DIMENSION(:  ), POINTER     :: min_phys,max_phys
    REAL(MK), DIMENSION(:  ), POINTER     :: cost

    INTEGER                               :: decomp,assig,tolexp
    INTEGER                               :: info,comm,rank,nproc
    INTEGER                               :: topoid
    INTEGER                               :: Npart=10000
    INTEGER                               :: i
    INTEGER,  DIMENSION(2*ndim)           :: bcdef
    INTEGER                               :: seedsize
    INTEGER,  DIMENSION(:),   ALLOCATABLE :: seed

    init
        USE ppm_module_init
#ifdef __MPI
        comm = MPI_COMM_WORLD
        CALL MPI_Comm_rank(comm,rank,info)
        CALL MPI_Comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = INT(LOG10(EPSILON(1.0_MK)))
        CALL ppm_init(ndim,MK,tolexp,0,debug,info,99)

        ALLOCATE(min_phys(ndim),max_phys(ndim),STAT=info)

        min_phys = 0.0_MK
        max_phys = 1.0_MK

        bcdef    = ppm_param_bcdef_freespace

        NULLIFY(xp)
    end init


    finalize
        USE ppm_module_finalize

        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys)
    end finalize


    setup

        CALL RANDOM_SEED(SIZE=seedsize)
        ALLOCATE(seed(seedsize),STAT=info)
        DO i=1,seedsize
           seed(i)=10+i*i*(rank+1)
        ENDDO
        CALL RANDOM_SEED(PUT=seed)

    end setup


    teardown

        DEALLOCATE(seed)

    end teardown

    test color_part
        ! test to check the vertex coloring results

        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_write
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_particles_typedef

        TYPE(ppm_t_particles_d), TARGET :: Part

        start_subroutine("color_part")

        CALL Part%create(Npart,info,name="Part1")
        Assert_Equal(info,0)

        CALL Part%get_xp(xp,info)
        Assert_Equal(info,0)

        !----------------
        ! create particles
        !----------------
        CALL RANDOM_NUMBER(xp)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_bisection
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,Npart,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
        Assert_Equal(info,0)

        CALL Part%set_xp(xp,info)
        Assert_Equal(info,0)

        Part%active_topoid=topoid
        Part%ghostlayer=cutoff

        CALL Part%apply_bc(info)
        Assert_Equal(info,0)

        CALL Part%map(info,global=.TRUE.,topoid=topoid)
        Assert_Equal(info,0)

        CALL Part%map_ghosts(info)
        Assert_Equal(info,0)

        CALL Part%destroy(info)
        Assert_Equal(info,0)

        ASSOCIATE (topo => ppm_topo(topoid)%t)
           stdout_f('(2(A,I0))',"procs color= ",'topo%ineighcolor(0)'," - max color=",'topo%ineighcolor(topo%nneighproc+1)')
        END ASSOCIATE

        end_subroutine()

    end test

end test_suite
