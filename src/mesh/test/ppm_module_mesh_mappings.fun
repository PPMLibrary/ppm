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
logical                         :: unifpatch

integer                         :: i,j,k,sizex,sizey
integer                         :: nsublist
integer, dimension(:  ),pointer :: isublist => NULL()
integer, dimension(2*ndim)      :: bcdef
integer                         :: bcdefX,bcdefY
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


!Uncomment for more thorough testing...
    !test ghost_mappings({decomp: [ppm_param_decomp_bisection,ppm_param_decomp_cuboid,ppm_param_decomp_xpencil, ppm_param_decomp_ypencil,ppm_param_decomp_xy_slab,ppm_param_decomp_bisection]}, sizex:[69,70,71], sizey:[77,78,81,82]})
    test ghost_mappings_basics({decomp: [ppm_param_decomp_cuboid]}, sizex:[69,70], sizey:[77,78,81,82],unifpatch: [.true.,.false.]})


        type(ppm_t_field) :: Field1,Field2
        real(ppm_kind_double),dimension(ndim) :: pos
        integer                             :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER    :: dinfo => NULL()
        logical                             :: assoc

        start_subroutine("ghost_mappings_basics")


        if (decomp.eq.ppm_param_decomp_xpencil .and. (sizey/nproc).LE.2) return
        if (decomp.eq.ppm_param_decomp_ypencil .and. (sizex/nproc).LE.2) return
        if (ndim.eq.2 .and. decomp.eq.ppm_param_decomp_xy_slab) return

        if (ppm_debug.ge.1 .and. rank.eq.0) then
            stdout("STARTING test with decomp = ",decomp,topoid,sizex,sizey)
        endif
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        offset = 0._mk
        Nm(ndim) = 65
        Nm(1:2) = (/sizex,sizey/)

        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = REAL(MAXVAL(ighostsize),MK)/MAXVAL(Nm)

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
            Assert_Equal(info,0)
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)

        IF (.NOT.unifpatch) THEN
            if (ndim.eq.2) then
                my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
                my_patch(1:4) = (/0.17_mk,0.17_mk,0.70_mk,0.70_mk/)
            else
                my_patch(1:6) = (/0.15_mk,0.10_mk,0.51_mk,&
                    0.99_mk,0.7_mk,0.78_mk/)
            endif
            call Mesh1%def_patch(my_patch,info)
            Assert_Equal(info,0)
        ENDIF

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
            Assert_True(ALL(p%istart.GE.0))
            Assert_True(ALL(p%iend.GE.0))
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
                DO ipatch=1,Mesh1%subpatch_by_sub(isub)%nsubpatch
                    SELECT TYPE(p => Mesh1%subpatch_by_sub(isub)%vec(ipatch)%t)
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
                DO ipatch=1,Mesh1%subpatch_by_sub(isub)%nsubpatch
                    SELECT TYPE(p => Mesh1%subpatch_by_sub(isub)%vec(ipatch)%t)
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
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
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
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        if (ppm_debug.ge.1 .and. rank.eq.0) then
            stdout("FINISHED test with decomp = ",decomp,topoid)
        endif
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        end_subroutine()
    end test

    !test ghost_mappings_bcdef({decomp: [ppm_param_decomp_cuboid], sizex:[69], sizey:[77], bcdefX: [ppm_param_bcdef_periodic,ppm_param_bcdef_symmetry,ppm_param_bcdef_freespace], bcdefY: [ppm_param_bcdef_periodic,ppm_param_bcdef_symmetry,ppm_param_bcdef_freespace],unifpatch: [.true.,.false.]})
    test ghost_mappings_bcdef({decomp: [ppm_param_decomp_cuboid], sizex:[69,70,71,72,73,75,123], sizey:[77,83,93,111], bcdefX: [ppm_param_bcdef_periodic,ppm_param_bcdef_freespace], bcdefY: [ppm_param_bcdef_periodic,ppm_param_bcdef_freespace],unifpatch: [.true.,.false.]})
        type(ppm_t_field) :: Field1,Field2
        real(ppm_kind_double),dimension(ndim) :: pos
        integer                             :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER    :: dinfo => NULL()
        logical                             :: assoc

        start_subroutine("ghost_mappings_bcdef")

        if (decomp.eq.ppm_param_decomp_xpencil .and. (sizey/nproc).LE.2) return
        if (decomp.eq.ppm_param_decomp_ypencil .and. (sizex/nproc).LE.2) return
        if (ndim.eq.2 .and. decomp.eq.ppm_param_decomp_xy_slab) return

        bcdef = ppm_param_bcdef_periodic
        bcdef(1) = bcdefX
        bcdef(2) = bcdefX
        bcdef(3) = bcdefY
        bcdef(4) = bcdefY

        if (ppm_debug.ge.1 .and. rank.eq.0) then
            stdout("STARTING test with bcdef = ",'bcdef(1:2*ndim)')
            stdout(" patch covers the whole domain? ",unifpatch)
        endif

        offset = 0._mk
        Nm(ndim) = 65
        Nm(1:2) = (/sizex,sizey/)


        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = REAL(MAXVAL(ighostsize),MK)/MAXVAL(Nm)

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
            Assert_Equal(info,0)
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)

        if (.not.unifpatch) then
            if (ndim.eq.2) then
                my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
                my_patch(1:4) = (/0.17_mk,0.17_mk,0.70_mk,0.70_mk/)
            else
                my_patch(1:6) = (/0.15_mk,0.10_mk,0.51_mk,0.99_mk,0.7_mk,0.78_mk/)
            endif

            call Mesh1%def_patch(my_patch,info)
            Assert_Equal(info,0)
        endif

        call Field1%create(2,info,name='vecField')
            Assert_Equal(info,0)
        call Field1%discretize_on(Mesh1,info)
            Assert_Equal(info,0)
        call Field2%create(1,info,name='scaField')
            Assert_Equal(info,0)
        call Field2%discretize_on(Mesh1,info)
            Assert_Equal(info,0)

        if (ppm_debug.GE.1) then
            topo => ppm_topo(Mesh1%topoid)%t
#ifdef __MPI
            call MPI_BARRIER(comm,info)
#endif
            stdout("NB subdomains =  ",topo%nsubs)
            do i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("********************************")
                stdout("Subdomain nb ",isub)
                stdout("coordinates subs Min =  ",'topo%min_subd(1:ndim,isub)')
                stdout("coordinates subs Max =  ",'topo%max_subd(1:ndim,isub)')
                stdout("mesh coords subs Min =  ",'Mesh1%istart(1:ndim,isub)')
                stdout("mesh coords subs Max =  ",'Mesh1%iend(1:ndim,isub)')
                stdout("********************************")
            enddo
#ifdef __MPI
            call MPI_BARRIER(comm,info)
#endif
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
                stdout("patch bc    are  ",'p%bc(1:2*ndim)')
                stdout("--------------------------------")
                p => Mesh1%subpatch%next()
            enddo
#ifdef __MPI
            call MPI_BARRIER(comm,info)
#endif
        endif


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
            Assert_True(ALL(p%istart.GE.0))
            Assert_True(ALL(p%iend.GE.0))
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
                DO ipatch=1,Mesh1%subpatch_by_sub(isub)%nsubpatch
                    SELECT TYPE(p => Mesh1%subpatch_by_sub(isub)%vec(ipatch)%t)
                    TYPE IS (ppm_t_subpatch)
                        DO j=1,p%nnodes(2)
                        DO i=1,p%nnodes(1)
                            pos(1:ndim) = p%get_pos(i,j)
                            Assert_True(ALL(pos(1:ndim).LE.topo%max_subd(1:ndim,isub)))
                            Assert_True(ALL(pos(1:ndim).GE.topo%min_subd(1:ndim,isub)))
                            IF (p%bc(2).EQ.ppm_param_bcdef_periodic) THEN
                                Assert_True(pos(1).LT.topo%max_subd(1,isub))
                            ENDIF
                            IF (p%bc(4).EQ.ppm_param_bcdef_periodic) THEN
                                Assert_True(pos(2).LT.topo%max_subd(2,isub))
                            ENDIF
                        ENDDO
                        ENDDO
                    END SELECT
                ENDDO
            ENDDO
        ELSE
            DO jsub = 1,topo%nsublist
                isub = topo%isublist(jsub)
                DO ipatch=1,Mesh1%subpatch_by_sub(isub)%nsubpatch
                    SELECT TYPE(p => Mesh1%subpatch_by_sub(isub)%vec(ipatch)%t)
                    TYPE IS (ppm_t_subpatch)
                        DO k=1,p%nnodes(3)
                        DO j=1,p%nnodes(2)
                        DO i=1,p%nnodes(1)
                            pos(1:ndim) = p%get_pos(i,j,k)
                            Assert_True(ALL(pos(1:ndim).LE.topo%max_subd(1:ndim,isub)))
                            Assert_True(ALL(pos(1:ndim).GE.topo%min_subd(1:ndim,isub)))
                            IF (p%bc(2).EQ.ppm_param_bcdef_periodic) THEN
                                Assert_True(pos(1).LT.topo%max_subd(1,isub))
                            ENDIF
                            IF (p%bc(4).EQ.ppm_param_bcdef_periodic) THEN
                                Assert_True(pos(2).LT.topo%max_subd(2,isub))
                            ENDIF
                            IF (p%bc(6).EQ.ppm_param_bcdef_periodic) THEN
                                Assert_True(pos(ndim).LT.topo%max_subd(3,isub))
                            ENDIF
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
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif

        call Mesh1%map_send(info)
            Assert_Equal(info,0)
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        call Field2%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)
        call Field1%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif

        !Now check that the ghost mapping has been done correctly
        ! by comparing the values of all nodes (incl. ghosts) to the
        ! theoretical values.
        nb_errors = 0
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for valid_nodes
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .lt. 0._mk) then
                    nb_errors = nb_errors + 1
                ENDIF
                !stdout("i,j = ",i,j,'sbpitr%istart(1)-1+i','sbpitr%istart(2)-1+j')
                !stdout("pos = ",pos)
                Assert_Equal_Within(Field1_n(1) ,cos(2._mk*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,cos(2._mk*pi*pos(1)) + 2._mk,1e-5)
                Assert_Equal_Within(Field2_n    ,1._mk,1e-5)
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for valid_nodes
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

#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif

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

#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        if (ppm_debug.ge.1 .and. rank.eq.0) then
            stdout("FINISHED test with decomp = ",decomp,topoid)
        endif
#ifdef __MPI
        call MPI_BARRIER(comm,info)
#endif
        end_subroutine()
    end test

    FUNCTION my_init_function(x) RESULT(val)
        REAL(ppm_kind_double) :: val
        REAL(ppm_kind_double),DIMENSION(1:ndim),INTENT(IN) :: x

        val = sum(x)
    END FUNCTION

end test_suite
