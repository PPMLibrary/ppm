test_suite ppm_module_mesh_small

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
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
integer                          :: ipatch,isub
class(ppm_t_subpatch_),POINTER   :: p => NULL()

integer                          :: mypatchid
real(mk),dimension(2*ndim)       :: my_patch
real(mk),dimension(ndim)         :: offset

real(mk),dimension(:,:),pointer  :: field2d_1,field2d_2
real(mk),dimension(:,:,:),pointer:: field3d_1,field3d_2
real(mk),dimension(:,:,:,:),pointer:: field4d_1,field4d_2

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
    test mesh_create_destroy
        start_subroutine("test_create_destroy")
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.05_mk 
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)




        call MPI_BARRIER(comm,info)
        topo => ppm_topo(topoid)%t
        stdout("NB subdomains GLOBAL =  ",topo%nsubs)
        stdout("NB subdomains LOCAL  =  ",topo%nsublist)
        do i = 1,topo%nsublist
            isub = topo%isublist(i)
            stdout("Subdomain number ",isub)
            stdout("coordinates  =  ",'topo%min_subd(1:2,isub)')
            stdout("number of neighs  =  ",'topo%nneighsubs(i)')
            do j=1,topo%nneighsubs(i)
            stdout("list of neighs        :  ",'topo%ineighsubs(j,i)')
            enddo
        enddo
        call MPI_BARRIER(comm,info)

        Nm = 125
        offset = 0._mk
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)
        call Mesh1%destroy(info)
        Assert_Equal(info,0)

        h = (max_phys-min_phys)/Nm
        call Mesh1%create(topoid,offset,info,h=h)
        Assert_Equal(info,0)
        call Mesh1%destroy(info)
        Assert_Equal(info,0)
        end_subroutine()
    end test


    test small_test
        use ppm_module_io_vtk
        type(ppm_t_field) :: Field1,Field2
        real(ppm_kind_double),dimension(ndim) :: pos
        real(mk),dimension(:,:,:),pointer :: Field1_data => NULL()
        real(mk),dimension(:,:),  pointer :: Field2_data => NULL()
        integer :: p_idx
        CLASS(ppm_t_discr_info_),POINTER             :: dinfo => NULL()
        logical :: assoc

        start_subroutine("test_patch_fields")

        Nm = 50
        Nm(ndim) = 50
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)

        my_patch(1:2*ndim) = (/0.15,0.10,0.8,0.7/)
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

        p_idx = Field1%get_pid(Mesh1)
        Assert_True(p_idx.GT.0)
        p_idx = Field2%get_pid(Mesh1)
        Assert_True(p_idx.GT.0)

        p => Mesh1%subpatch%begin()
        do while (ASSOCIATED(p))
            Assert_True(associated(p%mesh,Mesh1))
            Assert_True(associated(p%subpatch_data))
            p_idx = Field1%get_pid(p%mesh)
            Assert_True(p_idx.GT.0)
            p_idx = Field2%get_pid(p%mesh)
            Assert_True(p_idx.GT.0)
            p => Mesh1%subpatch%next()
        enddo

        if (ppm_debug.GT.0) then
            call MPI_BARRIER(comm,info)
            topo => ppm_topo(Mesh1%topoid)%t
            stdout("NB subdomains =  ",topo%nsubs)
            do i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("coordinates subs =  ",'topo%min_subd(1:2,isub)')
            enddo
            call MPI_BARRIER(comm,info)
            stdout("NB subpatch =  ",Mesh1%subpatch%nb)
            p => Mesh1%subpatch%begin()
            do while (ASSOCIATED(p))
                stdout("patch global ",'p%istart(1:2)',"/",'p%iend(1:2)')
                p => Mesh1%subpatch%next()
            enddo
            call MPI_BARRIER(comm,info)
        endif

        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                Field1_n(1) = -10.17_mk
                Field1_n(2) = -20.17_mk
                Field2_n    = -40.17_mk
        end foreach

        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
                Field1_n(1) = cos(2._mk*pi*pos(1))
                Field1_n(2) = cos(2._mk*pi*pos(1)) + 2._mk
                Field2_n    = cos(2._mk*pi*pos(1)) + 6._mk
        end foreach

        call Mesh1%map_ghost_get(info)
            Assert_Equal(info,0)

        call FLUSH()
        call MPI_BARRIER(comm,info)
        stdout("Ghost Get Done")
        call FLUSH()
        call MPI_BARRIER(comm,info)

        call Field1%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)

        call FLUSH()
        call MPI_BARRIER(comm,info)
        stdout("Ghost Push Done")
        call FLUSH()
        call MPI_BARRIER(comm,info)

        call Mesh1%map_send(info)
            Assert_Equal(info,0)

        call FLUSH()
        call MPI_BARRIER(comm,info)
        stdout("Send Done")
        call FLUSH()
        call MPI_BARRIER(comm,info)


        call Field1%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)

        call FLUSH()
        call MPI_BARRIER(comm,info)
        stdout("Ghost Pop Done")
        call FLUSH()
        call MPI_BARRIER(comm,info)



        call Mesh1%destroy(info)
            Assert_Equal(info,0)
        call Field1%destroy(info)
            Assert_Equal(info,0)
        call Field2%destroy(info)
            Assert_Equal(info,0)

        end_subroutine()
    end test


end test_suite
