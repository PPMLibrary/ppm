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
integer,parameter               :: ndim=3
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
integer                         :: topoid,nneigh_theo
integer                         :: np_global = 10000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.15_mk
real(mk)                        :: cutoff_input
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: rcp,wp
integer,dimension(:),pointer    :: wpi=>NULL()
integer                         :: i,j,k,isum1,isum2,ip,wp_id
integer                         :: wp1_id = 0, dwp1_id = 0
integer                         :: wp2_id = 0
integer                         :: wp3_id = 0
integer                         :: op_id
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
type(ppm_t_sop_d)               :: Pc_a
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()
real(mk)                        :: err

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


    test operators

        call Pc%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)

        call Pc%map(info,global=.true.,topoid=topoid)
        call Pc%map_ghosts(info)

        call Pc%comp_neighlist(info)
        Assert_Equal(info,0)

        call Pc%create_op(op_id,1,(/2._mk/),(/2,0,0/),(/2/),info,name='dx2')
        Assert_Equal(info,0)

        call Pc%comp_op(op_id,info)
        Assert_Equal(info,0)

        call Pc%create_prop(wp1_id,ppm_type_real_double,info,1,name='testf_sca')
        call Pc%create_prop(wp2_id,ppm_type_real_double,info,3,name='testf_vec')
        call Pc%get(wp_1r,wp1_id)
        call Pc%get(wp_2r,wp2_id)
        call Pc%get_xp(xp)
        DO ip=1,Pc%Npart
            wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
            wp_2r(1:ndim,ip) = f0_test(xp(1:ndim,ip),ndim)
        ENDDO
        call Pc%set_xp(xp,read_only=.true.)
        call Pc%set(wp_1r,wp1_id)
        call Pc%set(wp_2r,wp2_id)
        call Pc%map_ghosts(info)

        call Pc%apply_op(wp1_id,dwp1_id,op_id,info)
        Assert_Equal(info,0)

        call Pc%destroy_op(op_id,info)
        Assert_Equal(info,0)


    end test

    test initialize_cart
        ! test initialization of particles on a grid

        use ppm_module_typedef

        call Pc%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)


        call Pc%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Pc%map_ghosts(info)
        Assert_Equal(info,0)

        !creating/destroying properties of different types
        call Pc%create_prop(wp1_id,ppm_type_longint,info,&
            name='test_li',with_ghosts=.true.)
        Assert_Equal(info,0)
        call Pc%get(wp_1li,wp1_id,with_ghosts=.true.)
        call Pc%set(wp_1li,wp1_id)
        call Pc%destroy_prop(wp1_id,info)

        call Pc%create_prop(wp1_id,ppm_type_int,info,name='test_i',zero=.true.)
        Assert_Equal(info,0)
        call Pc%get(wp_1i,wp1_id)
        call Pc%set(wp_1i,wp1_id)

        call Pc%create_prop(wp1_id,ppm_type_int,info,10,&
            name='test_2i',zero=.true.)
        Assert_Equal(info,0)
        call Pc%get(wp_2i,wp1_id)
        Assert_Equal(MINVAL(wp_2i),0)
        Assert_Equal(MAXVAL(wp_2i),0)
        call Pc%set(wp_2i,wp1_id)


        call Pc%create_prop(wp1_id,ppm_type_real_double,info,&
            name='test_r',zero=.true.)
        Assert_Equal(info,0)
        call Pc%get(wp_1r,wp1_id)
        call Pc%set(wp_1r,wp1_id)

        call Pc%create_prop(wp1_id,ppm_type_logical,info,name='test_l',zero=.true.)
        Assert_Equal(info,0)
        call Pc%get(wp_1l,wp1_id)
        call Pc%set(wp_1l,wp1_id)


        DO i=1,25
            call Pc%create_prop(wp2_id,ppm_type_comp_double,info,3)
            Assert_Equal(info,0)
            call Pc%create_prop(wp3_id,ppm_type_real_double,info,1)
            Assert_Equal(info,0)
            call Pc%get(wp_2c,wp2_id)
            call Pc%set(wp_2c,wp2_id)

            call Pc%destroy_prop(wp3_id,info)
            Assert_Equal(info,0)
            call Pc%destroy_prop(wp2_id,info)
            Assert_Equal(info,0)
        ENDDO

        !Set up a velocity field and a scalar test function on the particles
        call Pc%create_prop(wp2_id,ppm_type_real_double,info,3,name='velocity')
        Assert_Equal(info,0)
        call Pc%create_prop(wp3_id,ppm_type_real_double,info,1,name='testf')
        Assert_Equal(info,0)
        call Pc%get(wp_2r,wp2_id)
        call Pc%get(wp_1r,wp3_id)
        call Pc%get_xp(xp)
        DO ip=1,Pc%Npart
            wp_2r(1:ndim,ip) = COS((10._MK*xp(1:ndim,ip))**2)
            wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
        ENDDO
        wp_2r = cos(wp_2r) * Pc%ghostlayer
        call Pc%set_xp(xp,read_only=.true.)
        call Pc%set(wp_1r,wp3_id)

        !Move the particles with this velocity field
        call Pc%move(wp_2r,info)
        Assert_Equal(info,0)

        call Pc%set(wp_2r,wp2_id)

        !Apply boundary conditions and remap the particles
        call Pc%apply_bc(info)
        Assert_Equal(info,0)
        call Pc%map(info)
        Assert_Equal(info,0)

        !Get the new ghosts
        call Pc%map_ghosts(info)
        Assert_Equal(info,0)


        !Move the particles back to their initial positions
        call Pc%get(wp_2r,wp2_id)
        wp_2r = -wp_2r
        call Pc%move(wp_2r,info)
        call Pc%set(wp_2r,wp2_id)


        !Re-apply boundary conditions and remap the particles
        call Pc%apply_bc(info)
        Assert_Equal(info,0)
        call Pc%map(info)
        Assert_Equal(info,0)

        call Pc%map_ghosts(info)
        Assert_Equal(info,0)

        !Compare values and check that they are still the same
        call Pc%get_xp(xp)
        call Pc%get(wp_1r,wp3_id)
        err = 0._mk
        DO ip=1,Pc%Npart
            err = max(err,abs(wp_1r(ip) - f0_test(xp(1:ndim,ip),ndim)))
        ENDDO
        Assert_Equal_Within(err,0,tol)
        call Pc%set_xp(xp,read_only=.true.)
        call Pc%set(wp_1r,wp3_id,read_only=.true.)



        !call Pc%print_info(info)
        !Assert_Equal(info,0)

    end test

    test neighlists
        ! test neighbour lists

        use ppm_module_typedef

        call Pc%destroy(info)
        call Pc%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)

        call Pc%map(info,global=.true.,topoid=topoid)
        call Pc%map_ghosts(info)

        call Pc%comp_neighlist(info)
        Assert_Equal(info,0)

        write(*,*) Pc%neighs%vec(1)%t%nneighmin
        write(*,*) Pc%neighs%vec(1)%t%nneighmax

    end test

    test sop_type
        ! test procedures for sop data structures

        use ppm_module_typedef

        write(*,*) 'SOP stuff:'
        call Pc_a%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)


        call Pc_a%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Pc_a%map_ghosts(info)
        Assert_Equal(info,0)

        call Pc_a%comp_neighlist(info)
        Assert_Equal(info,0)

        write(*,*) Pc_a%neighs%vec(1)%t%nneighmin
        write(*,*) Pc_a%neighs%vec(1)%t%nneighmax

        !creating/destroying properties of different types
        call Pc_a%create_prop(wp1_id,ppm_type_longint,info,&
            name='test_li',with_ghosts=.true.)
        Assert_Equal(info,0)
        call Pc_a%get(wp_1li,wp1_id,with_ghosts=.true.)
        call Pc_a%set(wp_1li,wp1_id)
        call Pc_a%destroy_prop(wp1_id,info)

        call Pc_a%create_prop(wp1_id,ppm_type_int,info,name='test_i',zero=.true.)
        Assert_Equal(info,0)
        call Pc_a%get(wp_1i,wp1_id)
        call Pc_a%set(wp_1i,wp1_id)

        !Set up a velocity field and a scalar test function on the particles
        call Pc_a%create_prop(wp2_id,ppm_type_real_double,info,3,name='velocity')
        Assert_Equal(info,0)
        call Pc_a%create_prop(wp3_id,ppm_type_real_double,info,1,name='testf')
        Assert_Equal(info,0)
        call Pc_a%get(wp_2r,wp2_id)
        call Pc_a%get(wp_1r,wp3_id)
        call Pc_a%get_xp(xp)
        DO ip=1,Pc_a%Npart
            wp_2r(1:ndim,ip) = COS((10._MK*xp(1:ndim,ip))**2)
            wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
        ENDDO
        wp_2r = cos(wp_2r) * Pc_a%ghostlayer
        call Pc_a%set_xp(xp,read_only=.true.)
        call Pc_a%set(wp_1r,wp3_id)

        !Move the particles with this velocity field
        call Pc_a%move(wp_2r,info)
        Assert_Equal(info,0)

        call Pc_a%set(wp_2r,wp2_id)

        !Apply boundary conditions and remap the particles
        call Pc_a%apply_bc(info)
        Assert_Equal(info,0)
        call Pc_a%map(info)
        Assert_Equal(info,0)

        !Get the new ghosts
        call Pc_a%map_ghosts(info)
        Assert_Equal(info,0)


        !Move the particles back to their initial positions
        call Pc_a%get(wp_2r,wp2_id)
        wp_2r = -wp_2r
        call Pc_a%move(wp_2r,info)
        call Pc_a%set(wp_2r,wp2_id)


        !Re-apply boundary conditions and remap the particles
        call Pc_a%apply_bc(info)
        Assert_Equal(info,0)
        call Pc_a%map(info)
        Assert_Equal(info,0)

        call Pc_a%map_ghosts(info)
        Assert_Equal(info,0)

        !Compare values and check that they are still the same
        call Pc_a%get_xp(xp)
        call Pc_a%get(wp_1r,wp3_id)
        err = 0._mk
        DO ip=1,Pc_a%Npart
            err = max(err,abs(wp_1r(ip) - f0_test(xp(1:ndim,ip),ndim)))
        ENDDO
        Assert_Equal_Within(err,0,tol)
        call Pc_a%set_xp(xp,read_only=.true.)
        call Pc_a%set(wp_1r,wp3_id,read_only=.true.)



        !call Pc_a%print_info(info)
        !Assert_Equal(info,0)

    end test
!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function f0_test(pos,ndim)

    real(mk)                              :: f0_test
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    f0_test =  sin(2._mk*pi*pos(1)) * cos(2._mk*pi*pos(2)) * &
        & sin(2._mk*pi*pos(ndim))

    return

end function f0_test

    

end test_suite
