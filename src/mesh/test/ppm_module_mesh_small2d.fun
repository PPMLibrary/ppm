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
integer                          :: ipatch,isub,jsub
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




        if (ppm_debug.gt.0) then
            call MPI_BARRIER(comm,info)
            topo => ppm_topo(topoid)%t
            stdout("NB subdomains GLOBAL =  ",topo%nsubs)
            stdout("NB subdomains LOCAL  =  ",topo%nsublist)
            do i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("Subdomain number ",isub)
                stdout("coordinates  Min =  ",'topo%min_subd(1:2,isub)')
                stdout("coordinates  Max =  ",'topo%max_subd(1:2,isub)')
                stdout("number of neighs  =  ",'topo%nneighsubs(i)')
                do j=1,topo%nneighsubs(i)
                stdout("list of neighs        :  ",'topo%ineighsubs(j,i)')
                enddo
            enddo
            call MPI_BARRIER(comm,info)
        endif

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
        integer :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER             :: dinfo => NULL()
        logical :: assoc

        start_subroutine("test_patch_fields")

        Nm = 25
        Nm(ndim) = 45
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)

        my_patch(1:2*ndim) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
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

        ! Check if the subpatch nodes are all exactly within the right subdomain
        topo => ppm_topo(Mesh1%topoid)%t
        ! loop through all subdomains on this processor
        DO jsub = 1,topo%nsublist
            isub = topo%isublist(jsub)
            DO ipatch=1,Mesh1%subpatch_by_sub(jsub)%nsubpatch
                SELECT TYPE(p => Mesh1%subpatch_by_sub(jsub)%vec(ipatch)%t)
                TYPE IS (ppm_t_subpatch)
                    DO j=1,p%nnodes(2)
                        DO i=1,p%nnodes(1)
                            pos(1:ndim) = p%get_pos(i,j)
                            Assert_True(ALL(pos(1:ndim).LT.topo%max_subd(1:ndim,isub)))
                            Assert_True(ALL(pos(1:ndim).GE.topo%min_subd(1:ndim,isub)))
                        ENDDO
                    ENDDO
                END SELECT
            ENDDO
        ENDDO

        if (ppm_debug.GT.0) then
            call MPI_BARRIER(comm,info)
            stdout("NB subdomains =  ",topo%nsubs)
            do i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("coordinates subs Min =  ",'topo%min_subd(1:2,isub)')
                stdout("coordinates subs Max =  ",'topo%max_subd(1:2,isub)')
            enddo
            call MPI_BARRIER(comm,info)
            stdout("NB patch =  ",Mesh1%npatch)
            stdout("NB subpatch =  ",Mesh1%subpatch%nb)
            p => Mesh1%subpatch%begin()
                if(associated(p)) then
                    stdout("********************************")
                    stdout("patch     istart_p ",'p%istart_p(1:2)')
                    stdout("patch     iend_p ",'p%iend_p(1:2)')
                    stdout("********************************")
                endif
            do while (ASSOCIATED(p))
                stdout("--------------------------------")
                stdout("patch     istart ",'p%istart(1:2)')
                stdout("patch     iend   ",'p%iend(1:2)')
                stdout("patch     g_start",'p%ghostsize(1)','p%ghostsize(3)')
                stdout("patch     g_end  ",'p%ghostsize(2)','p%ghostsize(4)')
                stdout("patch alloc from ",'p%lo_a(1:2)')
                stdout("            to   ",'p%hi_a(1:2)')
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
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .lt. 0._mk) then
                    nb_errors = nb_errors + 1
                ENDIF
                !stdout_f('(A,2(I0,1X),A,2(E24.16,1X))',&
                !    "Local  i,j ",i,j," Pos=",pos)
                !stdout_f('(A,2(I0,1X))',"Global i,j ",&
                !    'i+sbpitr%istart(1)-1','j+sbpitr%istart(2)-1')
                !stdout("p%istart=",'sbpitr%istart')
                !stdout("p%istart_pos=",'sbpitr%get_pos(1,1)')
                !stdout("p%iend=",'sbpitr%iend')
                !stdout("p%iend_pos=",&
                !    'sbpitr%get_pos(sbpitr%nnodes(1),sbpitr%nnodes(2))')
                Assert_Equal_Within(Field1_n(1) ,cos(2._mk*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,cos(2._mk*pi*pos(1)) + 2._mk,1e-5)
                Assert_Equal_Within(Field2_n    ,1._mk,1e-5)
        end foreach
        Assert_Equal(nb_errors,0)

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
