test_suite ppm_module_mesh

USE ppm_module_data
USE ppm_module_error
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_mesh_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_mktopo

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: pi = ACOS(-1.0_MK)
INTEGER,PARAMETER               :: ndim=3
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc
REAL(MK)                        :: tol
INTEGER                         :: topoid=-1
REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()

INTEGER, DIMENSION(:  ),POINTER :: ighostsize => NULL()
REAL(MK)                        :: sca_ghostsize

INTEGER                         :: i,j,k
INTEGER                         :: nsublist
INTEGER, DIMENSION(:  ),POINTER :: isublist => NULL()
INTEGER, DIMENSION(2*ndim)      :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
INTEGER, DIMENSION(:  ),POINTER :: nm => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: h => NULL()
TYPE(ppm_t_topo),       POINTER :: topo => NULL()

TYPE(ppm_t_equi_mesh),TARGET     :: Mesh1,Mesh2
INTEGER                          :: ipatch,isub,jsub
CLASS(ppm_t_subpatch_),POINTER   :: p => NULL()

INTEGER                          :: mypatchid
REAL(MK),DIMENSION(6)       :: my_patch
REAL(MK),DIMENSION(ndim)         :: offset

REAL(MK),DIMENSION(:,:),POINTER  :: field2d_1,field2d_2
REAL(MK),DIMENSION(:,:,:),POINTER:: field3d_1,field3d_2
REAL(MK),DIMENSION(:,:,:,:),POINTER:: field4d_1,field4d_2

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
        CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

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

    test mesh_create_destroy
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.05_MK
        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)

        Nm = 125
        offset = 0.0_MK
        CALL Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)
        CALL Mesh1%destroy(info)
        Assert_Equal(info,0)

        h = (max_phys-min_phys)/Nm
        CALL Mesh1%create(topoid,offset,info,h=h)
        Assert_Equal(info,0)
        CALL Mesh1%destroy(info)
        Assert_Equal(info,0)
    end test

    test mesh_add_patches
        INTEGER, DIMENSION(2) :: patchid
        Nm = 129
        Nm(ndim) = 29
        patchid = (/17,23/)

        CALL Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        !One patch
        my_patch(1:6) = (/0.5_MK,0.3_MK,0.1_MK,5.1_MK,1.1_MK,10.0_MK/)
        my_patch(ndim) = 5.1
        my_patch(2*ndim) = 10.0

        CALL Mesh1%def_patch(my_patch,info,patchid(1))
        Assert_Equal(info,0)
        Assert_True(ASSOCIATED(Mesh1%subpatch))

        ipatch = 0
        topo => ppm_topo(Mesh1%topoid)%t
        DO i = 1,topo%nsublist
            isub = topo%isublist(i)
            IF (all(my_patch(1:ndim).LT.topo%max_subd(1:ndim,isub)) .AND. &
                all(my_patch(ndim+1:2*ppm_dim).GT.topo%min_subd(1:ndim,isub)))&
                then
                !count one subpatch
                ipatch = ipatch + 1
            ENDIF
        ENDDO

        isub = 0
        p => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(p))
            isub = isub+1
            p => Mesh1%subpatch%next()
        ENDDO
        Assert_Equal(isub,ipatch)

        !Second patch
        my_patch(1:6) = (/0.5_MK,0.3_MK,0.1_MK,5.1_MK,1.1_MK,10.0_MK/)
        my_patch(ndim) = 5.1
        my_patch(2*ndim) = 10.0

        CALL Mesh1%def_patch(my_patch,info,patchid(2))
        Assert_Equal(info,0)


        !Check that the patchids have been set properly
        DO i=1,2
            Assert_Equal(Mesh1%patch%vec(i)%t%patchid,patchid(i))
        ENDDO

        CALL Mesh1%destroy(info)
        Assert_Equal(info,0)

        Assert_False(ASSOCIATED(Mesh1%subpatch))
    end test

    test mesh_add_many_patches
        start_subroutine("add_many_patches")
        Nm = 39
        Nm(ndim) = 129
        CALL Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        IF (ndim .EQ. 2) THEN

            mypatchid = 0
            DO i = 1,Nm(1)/4
            DO j = 1,Nm(2)/4
                mypatchid = mypatchid + 1
                my_patch(1:4) = (/h(1)*i,h(2)*j,h(1)*(i+4),h(2)*(j+4)/)
                CALL Mesh1%def_patch(my_patch,info,mypatchid)
                Assert_Equal(info,0)
                Assert_True(ASSOCIATED(Mesh1%subpatch))
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
                CALL Mesh1%def_patch(my_patch,info,mypatchid)
                Assert_Equal(info,0)
                Assert_True(ASSOCIATED(Mesh1%subpatch))
            ENDDO
            ENDDO
            ENDDO
        ENDIF

        CALL Mesh1%destroy(info)
        Assert_Equal(info,0)

        Assert_False(ASSOCIATED(Mesh1%subpatch))
        end_subroutine()
    end test

    test mesh_add_patch_uniform
        !testing for a single patch that covers the whole domain
        Nm = 129
        Nm(ndim) = 65
        CALL Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        topo => ppm_topo(Mesh1%topoid)%t

        mypatchid = 1
        my_patch(1:ndim)        = min_phys(1:ndim)
        my_patch(ndim+1:2*ndim) = max_phys(1:ndim)

        CALL Mesh1%def_patch(my_patch,info,mypatchid)
            Assert_Equal(info,0)

        Assert_True(ASSOCIATED(Mesh1%subpatch))

        isub = 0
        p => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(p))
            isub = isub+1
            p => Mesh1%subpatch%next()
        ENDDO
        Assert_Equal(isub,topo%nsublist)

        CALL Mesh1%destroy(info)
            Assert_Equal(info,0)

        !test the wrapper routine, which does the same thing
        CALL Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        CALL Mesh1%def_uniform(info)
            Assert_Equal(info,0)
        CALL Mesh1%destroy(info)
            Assert_Equal(info,0)
    end test

    test ghost_mappings
        IMPLICIT NONE

        TYPE(ppm_t_field) :: Field1,Field2
        REAL(ppm_kind_double),DIMENSION(ndim) :: pos
        INTEGER                             :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER    :: dinfo => NULL()
        LOGICAL                             :: assoc

        start_subroutine("test_ghost_mappings")

        Nm = 25
        Nm(ndim) = 45
        CALL Mesh1%create(topoid,offset,info,Nm=Nm,&
        & ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        IF (ndim.EQ.2) THEN
           my_patch(1:4) = (/0.15_MK,0.10_MK,0.99_MK,0.7_MK/)
        else
           my_patch(1:6) = (/0.15_MK,0.10_MK,0.41_MK,0.99_MK,0.7_MK,0.78_MK/)
        ENDIF

        CALL Mesh1%def_patch(my_patch,info)
        Assert_Equal(info,0)

        CALL Field1%create(2,info,name='vecField')
        Assert_Equal(info,0)
        CALL Field1%discretize_on(Mesh1,info)
        Assert_Equal(info,0)
        CALL Field2%create(1,info,name='scaField')
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
            Assert_True(ALL(p%istart_p.GE.0))
            Assert_True(ALL(p%iend_p.GE.0))
            p_idx = Field1%get_pid(p%mesh)
            Assert_True(p_idx.GT.0)
            p_idx = Field2%get_pid(p%mesh)
            Assert_True(p_idx.GT.0)
            p => Mesh1%subpatch%next()
        ENDDO

        ! Check if the subpatch nodes are all exactly within the right subdomain
        topo => ppm_topo(Mesh1%topoid)%t
        ! loop through all subdomains on this processor
        IF (ndim.EQ.2) THEN
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
                Field1_n(1) = -10.0_MK * (rank+1)
                Field1_n(2) = -10.0_MK * (rank+1) - 1
                Field2_n    = -1.0_MK
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for all
                Field1_n(1) = -10.0_MK * (rank+1)
                Field1_n(2) = -10.0_MK * (rank+1) - 1
                Field2_n    = -1.0_MK
        end foreach
        ENDIF

        !Change the values of the real nodes to something that depends
        ! on position (so that we can later check the values)
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
                Field1_n(1) = COS(2._MK*pi*pos(1))
                Field1_n(2) = COS(2._MK*pi*pos(1)) + 2.0_MK
                Field2_n    = 1.0_MK
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                Field1_n(1) = COS(2._MK*pi*pos(1))
                Field1_n(2) = COS(2._MK*pi*pos(1)) + 2.0_MK
                Field2_n    = 1.0_MK
        end foreach
        ENDIF

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

        !Now check that the ghost mapping has been done correctly
        ! by comparing the values of all nodes (incl. ghosts) to the
        ! theoretical values.
        nb_errors = 0
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .LT. 0.0_MK) THEN
                    nb_errors = nb_errors + 1
                ENDIF
                Assert_Equal_Within(Field1_n(1) ,COS(2._MK*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,COS(2._MK*pi*pos(1)) + 2.0_MK,1e-5)
                Assert_Equal_Within(Field2_n    ,1.0_MK,1e-5)
        end foreach
        ELSE
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
        ENDIF
        Assert_Equal(nb_errors,0)

        CALL Mesh1%destroy(info)
            Assert_Equal(info,0)
        CALL Field1%destroy(info)
            Assert_Equal(info,0)
        CALL Field2%destroy(info)
            Assert_Equal(info,0)

        end_subroutine()
    end test


    FUNCTION my_init_function(x) RESULT(val)
        REAL(ppm_kind_double) :: val
        REAL(ppm_kind_double),DIMENSION(1:ndim),INTENT(IN) :: x

        val = SUM(x)
    END FUNCTION
!============ Test cases ======================
!    test mesh_define
!        ! a simplistic test for checking if mesh_define is working
!
!        USE ppm_module_topo_typedef
!        USE ppm_module_mktopo
!        USE ppm_module_map_field
!        USE ppm_module_map_field_global
!        USE ppm_module_mesh_define
!        USE ppm_module_topo_get
!        USE ppm_module_topo_check
!
!
!        ALLOCATE(nm(ndim),STAT=info)
!        DO i=1,ndim
!            nm(i) = 32*nproc
!        ENDDO
!
!        !----------------
!        ! make topology
!        !----------------
!        decomp = ppm_param_decomp_cuboid
!        !decomp = ppm_param_decomp_xpencil
!        assig  = ppm_param_assign_internal
!
!        topoid = 0
!        meshid = -1
!
!        CALL ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,    &
!        &               bcdef,ghostsize,cost,nm,info)
!
!        CALL ppm_topo_get_meshinfo(topoid,meshid,nm,istart,ndata,maxndata,&
!                        isublist,nsublist,info)
!
!        ALLOCATE(field((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)),nsublist),&
!        &        STAT=info) ! 2d
!        ALLOCATE(field_ref((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)),nsublist),&
!        &        STAT=info) ! 2d
!
!!        ALLOCATE(field((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)), &
!!        &        (1-ghostsize(3)):(maxndata(3)+ghostsize(3)),nsublist),&
!!        &       STAT=info) ! 3d
!!        ALLOCATE(field_ref((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)), &
!!        &        (1-ghostsize(3)):(maxndata(3)+ghostsize(3)),nsublist),&
!!        &       STAT=info) ! 3d
!
!        DO i=1,ndim
!            h(i) = (max_phys(i) - min_phys(i)) / REAL(ndata(i,1)-1,mk)
!        ENDDO
!
!        meshid_ref = -1
!        nm_ref = nm
!        CALL ppm_mesh_define(topoid,meshid_ref,nm_ref,istart_ref,ndata_ref,info)
!        CALL ppm_map_field_global(topoid,topoid,meshid,meshid_ref,info)
!        CALL ppm_map_field_push(topoid,meshid,field,info)
!        CALL ppm_map_field_send(info)
!        CALL ppm_map_field_pop(topoid,meshid_ref,field_ref,ghostsize,info)
!
!
!        Assert_Equal(info,0)
!    end test


end test_suite
