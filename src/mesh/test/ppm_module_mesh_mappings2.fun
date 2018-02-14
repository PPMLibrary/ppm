test_suite ppm_module_mesh_mappings2

USE ppm_module_data
USE ppm_module_error
USE ppm_module_write
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_mesh_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_mktopo
USE ppm_module_topo_alloc

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: pi = ACOS(-1.0_MK)
INTEGER,PARAMETER               :: ndim=3
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc
REAL(MK)                        :: tol
INTEGER                         :: topoid=-1,topoid2=-1
REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()

INTEGER, DIMENSION(:  ),POINTER :: ighostsize => NULL()
REAL(MK)                        :: sca_ghostsize
LOGICAL                         :: unifpatch

INTEGER                         :: i,j,k,sizex,sizey
INTEGER                         :: nsublist
INTEGER, DIMENSION(:  ),POINTER :: isublist => NULL()
INTEGER, DIMENSION(2*ndim)      :: bcdef
INTEGER                         :: bcdefX,bcdefY
REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
INTEGER, DIMENSION(:  ),POINTER :: nm => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: h => NULL()
TYPE(ppm_t_topo),       POINTER :: topo => NULL()

TYPE(ppm_t_equi_mesh),TARGET     :: Mesh1,Mesh2
INTEGER                          :: ipatch,isub,jsub
CLASS(ppm_t_subpatch_),POINTER   :: p => NULL()

INTEGER                          :: mypatchid
REAL(MK),DIMENSION(6)            :: my_patch
REAL(MK),DIMENSION(ndim)         :: offset

!---------------- init -----------------------

    init

        USE ppm_module_topo_typedef
        USE ppm_module_init

        ALLOCATE(min_phys(ndim),max_phys(ndim),&
            &         ighostsize(ndim),nm(ndim),h(ndim))

        min_phys(1:ndim) = 0.0_MK
        max_phys(1:ndim) = 1.0_MK
        ighostsize(1:ndim) = 2
        bcdef(1:2*ndim) = ppm_param_bcdef_periodic
        tolexp = -12

#ifdef __MPI
        comm = MPI_COMM_WORLD
        CALL MPI_Comm_rank(comm,rank,info)
        CALL MPI_Comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        CALL ppm_init(ndim,mk,tolexp,comm,debug,info,99)

    end init

!----------------------------------------------

!---------------- finalzie --------------------


    finalize
        USE ppm_module_finalize

        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys,ighostsize,h,nm)

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


    test mappings

        CLASS(ppm_t_discr_info_),POINTER    :: dinfo => NULL()

        TYPE(ppm_t_field) :: Field1,Field2,Field3

        REAL(ppm_kind_double),DIMENSION(ndim) :: pos

        INTEGER :: p_idx, nb_errors

        LOGICAL :: assoc

        start_subroutine("ghost_mappings_basics")


        decomp=ppm_param_decomp_xy_slab

        IF (ppm_debug.GE.1 .AND. rank.EQ.0) THEN
           stdout("STARTING test with decomp = ",decomp,topoid,sizex,sizey)
        ENDIF
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        offset = 0.0_MK
        Nm = (/80,80,80/)

        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = REAL(MAXVAL(ighostsize),MK)/MAXVAL(Nm)

        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
        &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        CALL Mesh1%create(topoid,offset,info,Nm=Nm,&
        &    ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        CALL Mesh1%def_uniform(info)
        Assert_Equal(info,0)

        CALL Field1%create(2,info,name='vecField')
        Assert_Equal(info,0)
        CALL Field2%create(1,info,name='scaField')
        Assert_Equal(info,0)

        CALL Field1%discretize_on(Mesh1,info)
        Assert_Equal(info,0)
        CALL Field2%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        p => Mesh1%subpatch%begin()
        IF (ASSOCIATED(p)) THEN
           p_idx = Field1%get_pid(Mesh1)
           Assert_True(p_idx.GT.0)
           p_idx = Field2%get_pid(Mesh1)
           Assert_True(p_idx.GT.0)
        ENDIF
        DO WHILE (ASSOCIATED(p))
           Assert_True(ASSOCIATED(p%mesh,Mesh1))
           Assert_True(ASSOCIATED(p%subpatch_data))
           Assert_True(ALL(p%istart.GE.0))
           Assert_True(ALL(p%iend.GE.0))
           p_idx = Field1%get_pid(p%mesh)
           Assert_True(p_idx.GT.0)
           p_idx = Field2%get_pid(p%mesh)
           Assert_True(p_idx.GT.0)
           p => Mesh1%subpatch%next()
        ENDDO

        ! Check if the subpatch nodes are all exactly within the right subdomain
        topo => ppm_topo(Mesh1%topoid)%t
        ! loop through all subdomains on this processor
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

        !Fill in the allocated field arrays (incl. ghost nodes) with some data
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for all
                Field1_n(1) = -10.0_MK * (rank+1)
                Field1_n(2) = -10.0_MK * (rank+1) - 1
                Field2_n    = -1.0_MK
        end foreach

        !Change the values of the real nodes to something that depends
        ! on position (so that we can later check the values)
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                Field1_n(1) = COS(2._MK*pi*pos(1))
                Field1_n(2) = COS(2._MK*pi*pos(1)) + 2.0_MK
                Field2_n    = 1.0_MK
        end foreach

        !Do a ghost mapping
        CALL Mesh1%map_ghost_get(info)
        Assert_Equal(info,0)

        CALL Field1%map_ghost_push(Mesh1,info)
        Assert_Equal(info,0)
        CALL Field2%map_ghost_push(Mesh1,info)
        Assert_Equal(info,0)

        !CALL Mesh1%map_send(info)
        !non-blocking send
        CALL Mesh1%map_isend(info)
        Assert_Equal(info,0)

        CALL Field2%map_ghost_pop(Mesh1,info)
        Assert_Equal(info,0)
        CALL Field1%map_ghost_pop(Mesh1,info)
        Assert_Equal(info,0)
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        !Now check that the ghost mapping has been done correctly
        ! by comparing the values of all nodes (incl. ghosts) to the
        ! theoretical values.
        nb_errors = 0
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                IF (Field2_n .LT. 0.0_MK) THEN
                    nb_errors = nb_errors + 1
                ENDIF
                Assert_Equal_Within(Field1_n(1) ,COS(2._MK*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,COS(2._MK*pi*pos(1)) + 2.0_MK,1e-5)
                Assert_Equal_Within(Field2_n    ,1.0_MK,1e-5)
        end foreach
        Assert_Equal(nb_errors,0)

        !CALL Mesh1%print_vtk("Mesh1",info,Field=Field2)
        !Assert_Equal(info,0)

        sca_ghostsize = REAL(MAXVAL(ighostsize),MK)/MAXVAL(Nm)
        topoid2 = 0
        assig  = ppm_param_assign_internal
        decomp = ppm_param_decomp_cuboid

        CALL ppm_mktopo(topoid2,decomp,assig,min_phys,max_phys,&
        &    bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        CALL Mesh2%create(topoid2,offset,info,Nm=Nm,&
        &    ghostsize=ighostsize,name='Test_Mesh_2')
        Assert_Equal(info,0)

        CALL Mesh2%def_uniform(info)
        Assert_Equal(info,0)

        CALL Field2%discretize_on(Mesh2,info)
        Assert_Equal(info,0)

        !---------------------------------------------------------------------
        ! Global mapping to the new mesh in the new topology
        !---------------------------------------------------------------------
        CALL Mesh1%map(mesh2,info)
        Assert_Equal(info,0)

        CALL Field2%map_push(Mesh1,info)
        Assert_Equal(info,0)

        !CALL Mesh1%map_send(info)
        !non-blocking send
        CALL Mesh1%map_isend(info)
        Assert_Equal(info,0)

        CALL Field2%map_pop(Mesh2,info)
        Assert_Equal(info,0)

        CALL Mesh1%destroy(info)
        Assert_Equal(info,0)

        CALL Field1%destroy(info)
        Assert_Equal(info,0)

        ighostsize=1
        !for vtk output one layer of ghost will suffice

        !Do a ghost mapping for vtk output
        CALL Mesh2%map_ghost_get(info,ghostsize=ighostsize)
        Assert_Equal(info,0)

        CALL Field2%map_ghost_push(Mesh2,info)
        Assert_Equal(info,0)

        !CALL Mesh2%map_send(info)
        !non-blocking send
        CALL Mesh2%map_isend(info)
        Assert_Equal(info,0)

        CALL Field2%map_ghost_pop(Mesh2,info)
        Assert_Equal(info,0)
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
!        CALL Mesh2%print_vtk("Mesh2",info)
!        Assert_Equal(info,0)

        CALL Mesh2%destroy(info)
        Assert_Equal(info,0)

        CALL Field2%destroy(info)
        Assert_Equal(info,0)

        CALL ppm_topo_dealloc(ppm_topo(topoid)%t,info)
        Assert_Equal(info,0)
        CALL ppm_topo_dealloc(ppm_topo(topoid2)%t,info)
        Assert_Equal(info,0)
        DEALLOCATE(ppm_topo(topoid)%t,STAT=info)
        Assert_Equal(info,0)
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        IF (ppm_debug.GE.1 .AND. rank.EQ.0) THEN
            stdout("FINISHED test with decomp = ",decomp,topoid)
        ENDIF
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        end_subroutine()
    end test

    test selective_mappings
        USE ppm_module_util_commopt

        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        TYPE(ppm_t_field) :: Field1

        INTEGER, DIMENSION(:,:,:), POINTER :: wp
        INTEGER, DIMENSION(:), ALLOCATABLE :: sendlist,recvlist

        start_subroutine("selective_mappings")

        offset = 0.0_MK
        Nm = (/40,40,40/)
        ighostsize=1
        sca_ghostsize = MAXVAL(REAL(ighostsize,MK)/REAL(Nm-1,MK))

        decomp=ppm_param_decomp_xy_slab
        assig  = ppm_param_assign_internal
        topoid = 0

        ! create topology
        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys, &
        &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)

        ! communication optimization, which also take care of vertex coloring
        CALL ppm_util_commopt(topoid,info)
        Assert_Equal(info,0)

        topo => ppm_topo(topoid)%t

        ! create mesh
        CALL Mesh1%create(topoid,offset,info,Nm=Nm,&
        &    ghostsize=ighostsize,name='Mesh1')
        Assert_Equal(info,0)

        ! define uniform mesh everywhere
        CALL Mesh1%def_uniform(info)
        Assert_Equal(info,0)

        ! create scaler integer Field
        CALL Field1%create(1,info,dtype=ppm_type_int,name='scaField')
        Assert_Equal(info,0)

        ! descretize Field on the mesh
        CALL Field1%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        NULLIFY(wp)

        ! Initialize Field
        sbpitr => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(Field1,wp,info)
           Assert_Equal(info,0)
           IF (ppm_dim.EQ.3) THEN
              wp=99999
              wp(1,1,1)=ppm_rank
              wp(sbpitr%nnodes(1),1,1)=ppm_rank
              wp(1,sbpitr%nnodes(2),1)=ppm_rank
              wp(sbpitr%nnodes(1),sbpitr%nnodes(2),1)=ppm_rank
              wp(1,1,sbpitr%nnodes(3))=ppm_rank
              wp(sbpitr%nnodes(1),1,sbpitr%nnodes(3))=ppm_rank
              wp(1,sbpitr%nnodes(2),sbpitr%nnodes(3))=ppm_rank
              wp(sbpitr%nnodes(1),sbpitr%nnodes(2),sbpitr%nnodes(3))=ppm_rank
           ELSE
              wp=0
           ENDIF
           sbpitr => Mesh1%subpatch%next()
        ENDDO

        ! Define the selective sendlist and receivelist for selective ghost mappings
        ! This is done using processors coloring
        IF (topo%ineighcolor(0).EQ.1) THEN
           ALLOCATE(sendlist(COUNT(topo%ineighcolor(1:topo%nneighproc).NE.1)))
           ALLOCATE(recvlist(0))
           sendlist=PACK(topo%ineighproc(1:topo%nneighproc),topo%ineighcolor(1:topo%nneighproc).NE.1)
           recvlist=-1
        ELSE
           ALLOCATE(sendlist(0))
           ALLOCATE(recvlist(COUNT(topo%ineighcolor(1:topo%nneighproc).EQ.1)))
           sendlist=-1
           recvlist=PACK(topo%ineighproc(1:topo%nneighproc),topo%ineighcolor(1:topo%nneighproc).EQ.1)
        ENDIF

        ! ghost update
        CALL Mesh1%map_ghost_get(info,sendlist=sendlist,recvlist=recvlist)
        Assert_Equal(info,0)
        CALL Field1%map_ghost_push(Mesh1,info)
        Assert_Equal(info,0)
        CALL Mesh1%map_isend(info)
        Assert_Equal(info,0)
        CALL Field1%map_ghost_pop(Mesh1,info)
        Assert_Equal(info,0)

        ! Check the updated values
        sbpitr => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(Field1,wp,info)
           Assert_Equal(info,0)
           IF (ppm_dim.EQ.3) THEN
              IF (topo%ineighcolor(0).NE.1) THEN
                 Assert_True(wp(1,1,0).EQ.MOD(ppm_rank-1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(sbpitr%nnodes(1),1,0).EQ.MOD(ppm_rank-1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(1,sbpitr%nnodes(2),0).EQ.MOD(ppm_rank-1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(sbpitr%nnodes(1),sbpitr%nnodes(2),0).EQ.MOD(ppm_rank-1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(1,1,sbpitr%nnodes(3)+1).EQ.MOD(ppm_rank+1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(sbpitr%nnodes(1),1,sbpitr%nnodes(3)+1).EQ.MOD(ppm_rank+1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(1,sbpitr%nnodes(2),sbpitr%nnodes(3)+1).EQ.MOD(ppm_rank+1+ppm_nproc,ppm_nproc))
                 Assert_True(wp(sbpitr%nnodes(1),sbpitr%nnodes(2),sbpitr%nnodes(3)+1).EQ.MOD(ppm_rank+1+ppm_nproc,ppm_nproc))
              ENDIF
           ENDIF
           sbpitr => Mesh1%subpatch%next()
        ENDDO

        end_subroutine()
    end test

end test_suite
