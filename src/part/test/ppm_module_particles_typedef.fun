test_suite ppm_module_particles_typedef

use ppm_module_typedef
use ppm_module_data

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
integer                         :: topoid,nneigh_theo
integer                         :: np_global = 100000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.15_mk
real(mk)                        :: cutoff_input
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: rcp,wp
integer,dimension(:),pointer    :: wpi=>NULL()
integer                         :: i,j,k,isum1,isum2,ip,wp_id
integer                         :: wp1_id = 0
integer                         :: wp2_id = 0
integer                         :: wp3_id = 0
real(mk)                        :: rsum1,rsum2,rsum3,rsum4
integer                         :: nstep
real(mk),dimension(:),pointer   :: delta
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
character(len=ppm_char)         :: dirname
integer                         :: isymm = 0
logical                         :: lsymm = .false.,ok
real(mk)                        :: t0,t1,t2,t3
type(ppm_t_particles_d)         :: Pc
type(ppm_t_particles_d)         :: Pc2
type(ppm_t_particles_d)         :: Pc_cross
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()

integer, dimension(:), pointer                 :: wp_1i => NULL()
integer, dimension(:,:), pointer               :: wp_2i => NULL()
integer(ppm_kind_int64),dimension(:),  pointer :: wp_1li => NULL()
integer(ppm_kind_int64),dimension(:,:),pointer :: wp_2li => NULL()
real(mk), dimension(:),   pointer              :: wp_1r => NULL()
real(mk), dimension(:,:), pointer              :: wp_2r => NULL()
complex(mk), dimension(:),   pointer           :: wp_1c => NULL()
complex(mk), dimension(:,:), pointer           :: wp_2c => NULL()
logical, dimension(:),   pointer               :: wp_1l => NULL()

    init

        use ppm_module_typedef
        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         delta(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        
#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = int(log10(epsilon(1._mk)))+10
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys,delta)

    end finalize


    setup


    end setup
        

    teardown
        
        call Pc%destroy(info)

    end teardown

    test initialize_cart
        ! test initialization of particles on a grid

        use ppm_module_typedef

        call Pc%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)


        !creating/destroying properties of different types
        call Pc%create_prop(wp1_id,ppm_type_longint,info,name='test_li')
        Assert_Equal(info,0)
        call Pc%get(wp_1li,wp1_id)
        call Pc%set(wp_1li,wp1_id)
        call Pc%destroy_prop(wp1_id,info)

        wp1_id = 0
        call Pc%create_prop(wp1_id,ppm_type_int,info,name='test_i')
        Assert_Equal(info,0)
        call Pc%get(wp_1i,wp1_id)
        call Pc%set(wp_1i,wp1_id)

        wp1_id = 0
        call Pc%create_prop(wp1_id,ppm_type_int,info,10,&
            name='test_2i',zero=.true.)
        Assert_Equal(info,0)
        call Pc%get(wp_2i,wp1_id)
        Assert_Equal(MINVAL(wp_2i),0)
        Assert_Equal(MAXVAL(wp_2i),0)
        call Pc%set(wp_2i,wp1_id)

        wp1_id = 0
        call Pc%create_prop(wp1_id,ppm_type_real_double,info,name='test_r')
        Assert_Equal(info,0)
        call Pc%get(wp_1r,wp1_id)
        call Pc%set(wp_1r,wp1_id)

        wp1_id = 0
        call Pc%create_prop(wp1_id,ppm_type_logical,info,name='test_l')
        Assert_Equal(info,0)
        call Pc%get(wp_1l,wp1_id)
        call Pc%set(wp_1l,wp1_id)

        DO i=1,25
            wp2_id = 0
            call Pc%create_prop(wp2_id,ppm_type_comp_double,info,3)
            Assert_Equal(info,0)
            wp3_id = 0
            call Pc%create_prop(wp3_id,ppm_type_real_double,info,1)
            Assert_Equal(info,0)
            call Pc%destroy_prop(wp3_id,info)
            Assert_Equal(info,0)
        ENDDO

        call Pc%get(wp_2c,wp2_id)
        call Pc%set(wp_2c,wp2_id)

        call Pc%destroy(info)
        Assert_Equal(info,0)

    end test

end test_suite
