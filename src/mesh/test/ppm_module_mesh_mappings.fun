test_suite ppm_module_mesh_mappings

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_mktopo
use ppm_module_topo_alloc

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
real(mk),dimension(6)            :: my_patch
real(mk),dimension(ndim)         :: offset

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

    test ghost_map({decomp: [ppm_param_decomp_cuboid,ppm_param_decomp_xpencil, ppm_param_decomp_ypencil,ppm_param_decomp_xy_slab,ppm_param_decomp_bisection]})

        integer, dimension(2) :: patchid
        start_subroutine("ghost_map")

        if (ndim.eq.2 .and. decomp.eq.ppm_param_decomp_xy_slab) return

        Nm = 55
        Nm(ndim) = 55
        offset = 0._mk
        patchid = (/17,23/)

        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.05_mk 
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)


        h = (max_phys-min_phys)/Nm
        call Mesh1%create(topoid,offset,info,h=h)
        Assert_Equal(info,0)
        call Mesh1%destroy(info)
        Assert_Equal(info,0)

        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        !One patch
        my_patch(1:6) = (/0.5_mk,0.3_mk,0.1_mk,5.1_mk,1.1_mk,10.0_mk/)
        my_patch(ndim) = 5.1
        my_patch(2*ndim) = 10.0

        call Mesh1%def_patch(my_patch,info,patchid(1)) 
        Assert_Equal(info,0)
        Assert_True(associated(Mesh1%subpatch))

        ipatch = 0
        topo => ppm_topo(Mesh1%topoid)%t
        do i = 1,topo%nsublist
            isub = topo%isublist(i)
            if (all(my_patch(1:ndim).LT.topo%max_subd(1:ndim,isub)) .AND. &
                all(my_patch(ndim+1:2*ppm_dim).GT.topo%min_subd(1:ndim,isub)))&
                then
                !count one subpatch
                ipatch = ipatch + 1
            endif
        enddo   

        isub = 0
        p => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(p))
            isub = isub+1
            p => Mesh1%subpatch%next()
        ENDDO
        Assert_Equal(isub,ipatch)

        !Second patch
        my_patch(1:6) = (/0.5_mk,0.3_mk,0.1_mk,5.1_mk,1.1_mk,10.0_mk/)
        my_patch(ndim) = 5.1
        my_patch(2*ndim) = 10.0

        call Mesh1%def_patch(my_patch,info,patchid(2)) 
        Assert_Equal(info,0)

        
        !Check that the patchids have been set properly
        DO i=1,2
            Assert_Equal(Mesh1%patch%vec(i)%t%patchid,patchid(i))
        ENDDO


        if (ndim .eq. 2) then

            mypatchid = 0
            DO i = 1,Nm(1)/4
            DO j = 1,Nm(2)/4
                mypatchid = mypatchid + 1
                my_patch(1:4) = (/h(1)*i,h(2)*j,h(1)*(i+4),h(2)*(j+4)/)
                call Mesh1%def_patch(my_patch,info,mypatchid) 
                Assert_Equal(info,0)
                Assert_True(associated(Mesh1%subpatch))
            ENDDO
            ENDDO

        else
            mypatchid = 0
            DO i = 1,Nm(1)/4
            DO j = 1,Nm(2)/4
            DO k = 1,Nm(3)/4
                mypatchid = mypatchid + 1
                my_patch(1:6) = (/h(1)*i,h(2)*j,h(3)*k,&
                    h(1)*(i+4),h(2)*(j+4),h(3)*(k+4)/)
                call Mesh1%def_patch(my_patch,info,mypatchid) 
                Assert_Equal(info,0)
                Assert_True(associated(Mesh1%subpatch))
            ENDDO
            ENDDO
            ENDDO
        endif
        
        call Mesh1%destroy(info)
        Assert_Equal(info,0)

        Assert_False(associated(Mesh1%subpatch))

        call Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        topo => ppm_topo(Mesh1%topoid)%t

        mypatchid = 1
        my_patch(1:ndim)        = min_phys(1:ndim)
        my_patch(ndim+1:2*ndim) = max_phys(1:ndim)

        call Mesh1%def_patch(my_patch,info,mypatchid) 
            Assert_Equal(info,0)

        Assert_True(associated(Mesh1%subpatch))

        isub = 0
        p => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(p))
            isub = isub+1
            p => Mesh1%subpatch%next()
        ENDDO
        Assert_Equal(isub,topo%nsublist)

        call Mesh1%destroy(info)
            Assert_Equal(info,0)

        !test the wrapper routine, which does the same thing
        call Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        call Mesh1%def_uniform(info)
            Assert_Equal(info,0)
        call Mesh1%destroy(info)
            Assert_Equal(info,0)

        call ppm_topo_dealloc(topo,info)
            Assert_Equal(info,0)
    end test

    test ghost_mappings({decomp: [ppm_param_decomp_bisection,ppm_param_decomp_bisection,ppm_param_decomp_cuboid,ppm_param_decomp_xpencil, ppm_param_decomp_ypencil,ppm_param_decomp_xy_slab,ppm_param_decomp_bisection]})


        type(ppm_t_field) :: Field1,Field2
        real(ppm_kind_double),dimension(ndim) :: pos
        integer                             :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER    :: dinfo => NULL()
        logical                             :: assoc
        start_subroutine("ghost_mappings")

        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.05_mk 
        if (ndim.eq.2 .and. decomp.eq.ppm_param_decomp_xy_slab) return
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)

        offset = 0._mk
        Nm = 25
        Nm(ndim) = 45
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
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

         stdout("OK decomp = ",decomp,topoid)

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
        IF (ndim.eq.2) THEN
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
        ELSE
            DO jsub = 1,topo%nsublist
                isub = topo%isublist(jsub)
                DO ipatch=1,Mesh1%subpatch_by_sub(jsub)%nsubpatch
                    SELECT TYPE(p => Mesh1%subpatch_by_sub(jsub)%vec(ipatch)%t)
                    TYPE IS (ppm_t_subpatch)
                        DO k=1,p%nnodes(3)
                        DO j=1,p%nnodes(2)
                        DO i=1,p%nnodes(1)
                            pos(1:ndim) = p%get_pos(i,j,k)
                            Assert_True(ALL(pos(1:ndim).LT.topo%max_subd(1:ndim,isub)))
                            Assert_True(ALL(pos(1:ndim).GE.topo%min_subd(1:ndim,isub)))
                        ENDDO
                        ENDDO
                        ENDDO
                    END SELECT
                ENDDO
            ENDDO
        ENDIF

        !Fill in the allocated field arrays (incl. ghost nodes) with some data
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                Field1_n(1) = -10._mk * (rank+1)
                Field1_n(2) = -10._mk * (rank+1) - 1
                Field2_n    = -1._mk
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for all
                Field1_n(1) = -10._mk * (rank+1)
                Field1_n(2) = -10._mk * (rank+1) - 1
                Field2_n    = -1._mk
        end foreach
        ENDIF

        !Change the values of the real nodes to something that depends
        ! on position (so that we can later check the values)
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
                Field1_n(1) = cos(2._mk*pi*pos(1))
                Field1_n(2) = cos(2._mk*pi*pos(1)) + 2._mk
                Field2_n    = 1._mk
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                Field1_n(1) = cos(2._mk*pi*pos(1))
                Field1_n(2) = cos(2._mk*pi*pos(1)) + 2._mk
                Field2_n    = 1._mk
        end foreach
        ENDIF

        !Do a ghost mapping
        call Mesh1%map_ghost_get(info)
            Assert_Equal(info,0)

        call Field1%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)
        call Field2%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)

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
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .lt. 0._mk) then
                    nb_errors = nb_errors + 1
                ENDIF
                Assert_Equal_Within(Field1_n(1) ,cos(2._mk*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,cos(2._mk*pi*pos(1)) + 2._mk,1e-5)
                Assert_Equal_Within(Field2_n    ,1._mk,1e-5)
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                IF (Field2_n .lt. 0._mk) then
                    nb_errors = nb_errors + 1
                ENDIF
                Assert_Equal_Within(Field1_n(1) ,cos(2._mk*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,cos(2._mk*pi*pos(1)) + 2._mk,1e-5)
                Assert_Equal_Within(Field2_n    ,1._mk,1e-5)
        end foreach
        ENDIF
        Assert_Equal(nb_errors,0)

        call Mesh1%destroy(info)
            Assert_Equal(info,0)
        call Field1%destroy(info)
            Assert_Equal(info,0)
        call Field2%destroy(info)
            Assert_Equal(info,0)

        call ppm_topo_dealloc(ppm_topo(topoid)%t,info)
            Assert_Equal(info,0)
        deallocate(ppm_topo(topoid)%t,STAT=info)
            Assert_Equal(info,0)

        end_subroutine()
    end test


    FUNCTION my_init_function(x) RESULT(val)
        REAL(ppm_kind_double) :: val
        REAL(ppm_kind_double),DIMENSION(1:ndim),INTENT(IN) :: x

        val = sum(x)
    END FUNCTION

end test_suite
