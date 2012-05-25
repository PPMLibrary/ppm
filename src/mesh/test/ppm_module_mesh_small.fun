test_suite ppm_module_mesh_small

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_particles_typedef
use ppm_module_mktopo

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = ACOS(-1._mk)
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
real(mk)                        :: tol
integer                         :: topoid=-1
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()

integer, dimension(:  ),pointer :: ighostsize => NULL()
real(mk)                        :: sca_ghostsize

integer                         :: i,j,k
integer                         :: nsublist
integer, dimension(:  ),pointer :: isublist => NULL()
integer, dimension(2*ndim)      :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
type(ppm_t_topo),       pointer :: topo => NULL()

type(ppm_t_equi_mesh),TARGET     :: Mesh1,Mesh2
integer                          :: ipatch,isub,jsub
class(ppm_t_subpatch_),POINTER   :: p => NULL()

integer                          :: mypatchid
real(mk),dimension(2*ndim)       :: my_patch
real(mk),dimension(ndim)         :: offset

integer, dimension(:), pointer                 :: wp_1i => NULL()
integer, dimension(:,:), pointer               :: wp_2i => NULL()
integer(ppm_kind_int64),dimension(:),  pointer :: wp_1li => NULL()
integer(ppm_kind_int64),dimension(:,:),pointer :: wp_2li => NULL()
real(mk), dimension(:),   pointer              :: wp_1r => NULL()
real(mk), dimension(:,:), pointer              :: wp_2r => NULL()
complex(mk), dimension(:),   pointer           :: wp_1c => NULL()
complex(mk), dimension(:,:), pointer           :: wp_2c => NULL()
logical, dimension(:),   pointer               :: wp_1l => NULL()

!---------------- init -----------------------

    init

        use ppm_module_topo_typedef
        use ppm_module_init
        
        allocate(min_phys(ndim),max_phys(ndim),&
            &         ighostsize(ndim),nm(ndim),h(ndim))
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        ighostsize(1:ndim) = 2
        bcdef(1:2*ndim) = ppm_param_bcdef_periodic
        tolexp = -12

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

!----------------------------------------------

!---------------- finalzie --------------------


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,ighostsize,h,nm)

    end finalize

!----------------------------------------------

!------------- setup --------------------------

    setup

    end setup
!----------------------------------------------
        

!--------------- teardown ---------------------
    teardown
        NULLIFY(topo)
    end teardown
!----------------------------------------------

    test small_test
        type(ppm_t_field) :: Field1,Field2
        type(ppm_t_particles_d) :: Part1
        real(ppm_kind_double),dimension(ndim) :: pos
        integer :: p_idx, nb_errors
        integer :: np_global = 3000
        CLASS(ppm_t_discr_info_),POINTER              :: dinfo => NULL()
        real(ppm_kind_double),dimension(:  ), POINTER :: up_1d => NULL()
        real(ppm_kind_double),dimension(:,:), POINTER :: up_2d => NULL()
        logical :: assoc

        start_subroutine("test_interp")

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.07_mk 
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)


        Nm = 35
        Nm(ndim) = 45
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)

        call Part1%initialize(np_global,info,topoid=topoid,name="Part1")
            Assert_Equal(info,0)

        call Part1%set_cutoff(3._mk * Part1%h_avg,info)
        Assert_Equal(info,0)

        allocate(wp_2r(ndim,Part1%Npart))
        call random_number(wp_2r)
        wp_2r = (wp_2r - 0.5_mk) * Part1%h_avg * 0.15_mk
        call Part1%move(wp_2r,info)
        Assert_Equal(info,0)
        deallocate(wp_2r)

        Assert_true(Part1%has_neighlist(Part1))
        call Part1%apply_bc(info)
        Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)


        if (ndim.eq.2) then
            my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.51_mk,0.99_mk,0.7_mk,0.78_mk/)
        endif
        call Mesh1%def_patch(my_patch,info) 
        Assert_Equal(info,0)

        call Field1%create(2,info,name='vecField') 
            Assert_Equal(info,0)
        call Field1%discretize_on(Mesh1,info)
            Assert_Equal(info,0)
        call Field2%create(1,info,name='scaField') 
            Assert_Equal(info,0)
        call Field2%discretize_on(Mesh1,info)
            Assert_Equal(info,0)


        p => Mesh1%subpatch%begin()
        IF (associated(p)) THEN
            p_idx = Field1%get_pid(Mesh1)
            Assert_True(p_idx.GT.0)
            p_idx = Field2%get_pid(Mesh1)
            Assert_True(p_idx.GT.0)
        ENDIF
        do while (ASSOCIATED(p))
            Assert_True(associated(p%mesh,Mesh1))
            Assert_True(associated(p%subpatch_data))
            Assert_True(ALL(p%istart_p.GE.0))
            Assert_True(ALL(p%iend_p.GE.0))
            p_idx = Field1%get_pid(p%mesh)
            Assert_True(p_idx.GT.0)
            p_idx = Field2%get_pid(p%mesh)
            Assert_True(p_idx.GT.0)
            p => Mesh1%subpatch%next()
        enddo


        if (ppm_debug.GT.0) then
            call MPI_BARRIER(comm,info)
            stdout("NB subdomains =  ",topo%nsubs)
            do i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("coordinates subs Min =  ",'topo%min_subd(1:ndim,isub)')
                stdout("coordinates subs Max =  ",'topo%max_subd(1:ndim,isub)')
            enddo
            call MPI_BARRIER(comm,info)
            stdout("NB patch =  ",Mesh1%npatch)
            stdout("NB subpatch =  ",Mesh1%subpatch%nb)
            p => Mesh1%subpatch%begin()
                if(associated(p)) then
                    stdout("********************************")
                    stdout("patch     istart_p ",'p%istart_p(1:ndim)')
                    stdout("patch     iend_p ",'p%iend_p(1:ndim)')
                    stdout("********************************")
                endif
            do while (ASSOCIATED(p))
                stdout("--------------------------------")
                stdout("patch     istart ",'p%istart(1:ndim)')
                stdout("patch     iend   ",'p%iend(1:ndim)')
                stdout("patch alloc from ",'p%lo_a(1:ndim)')
                stdout("            to   ",'p%hi_a(1:ndim)')
                stdout("--------------------------------")
                p => Mesh1%subpatch%next()
            enddo
            call MPI_BARRIER(comm,info)
        endif

        !Fill in the allocated field arrays (incl. ghost nodes) with some data
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                Field1_n(1) = -10._mk * (rank+1)
                Field1_n(2) = -10._mk * (rank+1) - 1
                Field2_n    = -1._mk
        end foreach

        !Change the values of the real nodes to something that depends
        ! on position (so that we can later check the values)
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
                Field1_n(1) = cos(2._mk*pi*pos(1))
                Field1_n(2) = cos(2._mk*pi*pos(1)) + 2._mk
                Field2_n    = 1._mk
        end foreach

        !Do a ghost mapping
        call Mesh1%map_ghost_get(info)
            Assert_Equal(info,0)

        call Field1%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)
        call Field2%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)

        call MPI_BARRIER(comm,info)

        call Mesh1%map_send(info)
            Assert_Equal(info,0)

        call Field2%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)
        call Field1%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)

        !Now check that the ghost mapping has been done correctly
        ! by comparing the values of all nodes (incl. ghosts) to the
        ! theoretical values.
        nb_errors = 0
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real_and_ghosts
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .lt. 0._mk) then
                    nb_errors = nb_errors + 1
                ENDIF
                Assert_Equal_Within(Field1_n(1) ,cos(2._mk*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,cos(2._mk*pi*pos(1)) + 2._mk,1e-5)
                Assert_Equal_Within(Field2_n    ,1._mk,1e-5)
        end foreach
        Assert_Equal(nb_errors,0)

        call MPI_BARRIER(comm,info)

        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                Field1_n(1) = 2._mk 
                Field1_n(2) = 4._mk
                Field2_n    = -1._mk
        end foreach

        call Mesh1%interp_to_part(Part1,Field1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        !call Part1%map_ghosts(info)
        !    Assert_Equal(info,0)

        call Part1%get_field(Field1,up_2d,info)
            Assert_Equal(info,0)

            stdout("Min/Max 1 = ",'MINVAL(up_2d(1,1:Part1%Npart))',&
                'MAXVAL(up_2d(1,1:Part1%Npart))')
            stdout("Min/Max 2 = ",'MinVAL(up_2d(2,1:Part1%Npart))',&
                'MAXVAL(up_2d(2,1:Part1%Npart))')


        call Mesh1%destroy(info)
            Assert_Equal(info,0)
        call Field1%destroy(info)
            Assert_Equal(info,0)
        call Field2%destroy(info)
            Assert_Equal(info,0)
        call Part1%destroy(info)
            Assert_Equal(info,0)

        end_subroutine()
    end test


end test_suite
