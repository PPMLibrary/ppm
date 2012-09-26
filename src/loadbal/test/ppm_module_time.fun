test_suite ppm_module_time

use ppm_module_particles_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_operator_typedef
use ppm_module_sop
use ppm_module_interfaces
use ppm_module_data
use ppm_module_util_dbg
use ppm_module_time
use ppm_module_loadbal

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = ACOS(-1._mk)
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc,topoid
integer                         :: np_global = 3000
real(mk),parameter              :: cutoff = 0.15_mk
real(mk),dimension(:,:),pointer :: xp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
integer                         :: i,j,k,ip,wp_id
integer                         :: nstep
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
integer                         :: isymm = 0
real(mk)                        :: t0,t1,t2,t3
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

integer, dimension(:),allocatable              :: degree,order
real(ppm_kind_double),dimension(:),allocatable :: coeffs
integer                                        :: nterms

    init

        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
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

        deallocate(min_phys,max_phys,len_phys)

    end finalize


    setup


    end setup
        

    teardown
        
    end teardown

!    test initialize_random
!        ! test initialization of particles sampled unif. random
!        type(ppm_t_particles_d)               :: Part1
!
!        call Part1%initialize(np_global,info,topoid=topoid,&
!            distrib=ppm_param_part_init_random)
!        Assert_Equal(info,0)
!        Assert_True(associated(Part1%xp))
!
!        call Part1%get_xp(xp,info)
!        Assert_Equal(info,0)
!        call Part1%set_xp(xp,info,read_only=.true.)
!        Assert_Equal(info,0)
!
!        call Part1%map(info,global=.true.,topoid=topoid)
!        Assert_Equal(info,0)
!        call Part1%map_ghosts(info)
!        Assert_Equal(info,0)
!
!        call Part1%destroy(info)
!        Assert_Equal(info,0)
!
!    end test

    test initialize_cart
        ! test initialization of particles on a grid
        type(ppm_t_particles_d)               :: Part1
        class(ppm_t_discr_data),POINTER       :: Prop1=>NULL(),Prop2=>NULL()
        type(ppm_t_field)                     :: Field1
        INTEGER,DIMENSION(:),POINTER          :: nvlist => NULL()
        INTEGER,DIMENSION(:,:),POINTER        :: vlist => NULL()
        class(ppm_t_neighlist_d_),POINTER     :: Nlist => NULL()
        integer                               :: np_local
        REAL(mk), DIMENSION(:  ), POINTER       :: cost
        !!! Aggregate cost for each subdomain
        REAL(mk), DIMENSION(:  ), POINTER       :: proc_cost


        call Part1%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)
        Assert_True(associated(Part1%xp))

        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
        

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)
        print*,size(xp)/2
!        call Part1%apply_bc(info)
!        Assert_Equal(info,0)
        np_local = INT(size(xp)/2)
!        call ppm_dbg_print(topoid,0._mk,1,1,info,xp,np_local)
!        Assert_Equal(info,0)

        !!! CREATE NLIST AND ESTIMATE THE COST
        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)
        print*,'DEBUG =0='        
        Nlist => Part1%get_neighlist(Part1)
        Assert_true(associated(Nlist))
        
        Assert_True(Part1%has_neighlist())
        call Part1%get_vlist(nvlist,vlist,info,Nlist)
        Assert_true(associated(vlist))
        Assert_true(associated(nvlist))
        Assert_Equal(info,0)
!        print*,'DEBUG =3='
        print*,rank,'=> size of nvlist',size(nvlist)
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        call ppm_get_cost(topoid,-1,xp,Part1%Npart,cost,proc_cost,info,vlist,nvlist)
        Assert_Equal(info,0)
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
        
        !Set up a velocity field and a scalar test function on the particles
        call Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,&
            lda=3,name="velocity")
        Assert_Equal(info,0)
        call Part1%create_prop(info,discr_data=Prop2,dtype=ppm_type_real,&
            lda=1,name="testf")
        Assert_Equal(info,0)
        call Part1%get(Prop1,wp_2r,info)
        Assert_Equal(info,0)
        call Part1%get(Prop2,wp_1r,info)
        Assert_Equal(info,0)
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        DO ip=1,Part1%Npart
            wp_2r(1:ndim,ip) = COS((10._MK*xp(1:ndim,ip))**2)
            wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
        ENDDO
        wp_2r = cos(wp_2r) * Part1%ghostlayer

        !Move the particles with this displacement field
        call Part1%move(wp_2r,info)
        Assert_Equal(info,0)

        !Apply boundary conditions and remap the particles
        call Part1%apply_bc(info)
        Assert_Equal(info,0)
        call Part1%map(info)
        Assert_Equal(info,0)

        !Get the new ghosts
        call Part1%map_ghosts(info)
        Assert_Equal(info,0)
        
        !!! CREATE NLIST AND ESTIMATE THE COST
        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)
        print*,'DEBUG =0='        
        Nlist => Part1%get_neighlist(Part1)
        Assert_true(associated(Nlist))
        
        Assert_True(Part1%has_neighlist())
        call Part1%get_vlist(nvlist,vlist,info,Nlist)
        Assert_true(associated(vlist))
        Assert_true(associated(nvlist))
        Assert_Equal(info,0)
!        print*,'DEBUG =3='
        print*,rank,'=> size of nvlist',size(nvlist)
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        call ppm_get_cost(topoid,-1,xp,Part1%Npart,cost,proc_cost,info,vlist,nvlist)
        Assert_Equal(info,0)
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
        
        !!! MOVE THE PARTICLES BACK AND CHECK THE COSTS AGAIN
!
        !Move the particles back to their initial positions
        call Part1%get(Prop1,wp_2r,info)
        wp_2r = -wp_2r
        call Part1%move(wp_2r,info)
        Assert_Equal(info,0)

        !Re-apply boundary conditions and remap the particles
        call Part1%apply_bc(info)
        Assert_Equal(info,0)
        call Part1%map(info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)
        
        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)

        Assert_True(Part1%has_neighlist())

        !Compare values and check that they are still the same
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        call Part1%get(Prop2,wp_1r,info)
        Assert_Equal(info,0)
        err = 0._mk
        DO ip=1,Part1%Npart
            err = max(err,abs(wp_1r(ip) - f0_test(xp(1:ndim,ip),ndim)))
        ENDDO
        Assert_Equal_Within(err,0,tol)
!

!!! CREATE NLIST AND ESTIMATE THE COST
        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)
        print*,'DEBUG =0='        
        Nlist => Part1%get_neighlist(Part1)
        Assert_true(associated(Nlist))
        
        Assert_True(Part1%has_neighlist())
        call Part1%get_vlist(nvlist,vlist,info,Nlist)
        Assert_true(associated(vlist))
        Assert_true(associated(nvlist))
        Assert_Equal(info,0)
!        print*,'DEBUG =3='
        print*,rank,'=> size of nvlist',size(nvlist)
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        call ppm_get_cost(topoid,-1,xp,Part1%Npart,cost,proc_cost,info,vlist,nvlist)
        Assert_Equal(info,0)
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
!        !call Part1%print_info(info)
!        !Assert_Equal(info,0)
!        call Part1%destroy(info)
!        Assert_Equal(info,0)
!        call Field1%destroy(info)
!        Assert_Equal(info,0)
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
