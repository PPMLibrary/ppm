test_suite ppm_module_mesh_small2d

USE ppm_module_write
USE ppm_module_mesh_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_mktopo

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: pi = ACOS(-1.0_MK)
INTEGER,PARAMETER               :: ndim=2
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
        start_subroutine("test_create_destroy")
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




        IF (ppm_debug.GT.0) THEN
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
            topo => ppm_topo(topoid)%t
            stdout("NB subdomains GLOBAL =  ",topo%nsubs)
            stdout("NB subdomains LOCAL  =  ",topo%nsublist)
            DO i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("Subdomain number ",isub)
                stdout("coordinates  Min =  ",'topo%min_subd(1:2,isub)')
                stdout("coordinates  Max =  ",'topo%max_subd(1:2,isub)')
                stdout("number of neighs  =  ",'topo%nneighsubs(i)')
                DO j=1,topo%nneighsubs(i)
                stdout("list of neighs        :  ",'topo%ineighsubs(j,i)')
                ENDDO
            ENDDO
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
        ENDIF

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
        end_subroutine()
    end test

    test small_test
        TYPE(ppm_t_field) :: Field1,Field2
        REAL(ppm_kind_double),DIMENSION(ndim) :: pos
        REAL(MK),DIMENSION(:,:,:),POINTER :: Field1_data => NULL()
        REAL(MK),DIMENSION(:,:),  POINTER :: Field2_data => NULL()
        INTEGER :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER             :: dinfo => NULL()
        LOGICAL :: assoc

        start_subroutine("test_patch_fields")

        Nm = 25
        Nm(ndim) = 45
        CALL Mesh1%create(topoid,offset,info,Nm=Nm,&
        &    ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        my_patch(1:2*ndim) = (/0.15_MK,0.10_MK,0.99_MK,0.7_MK/)
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

        IF (ppm_debug.GT.0) THEN
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
            stdout("NB subdomains =  ",topo%nsubs)
            DO i = 1,topo%nsublist
                isub = topo%isublist(i)
                stdout("coordinates subs Min =  ",'topo%min_subd(1:2,isub)')
                stdout("coordinates subs Max =  ",'topo%max_subd(1:2,isub)')
            ENDDO
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
            stdout("NB patch =  ",Mesh1%npatch)
            stdout("NB subpatch =  ",Mesh1%subpatch%nb)
            p => Mesh1%subpatch%begin()
                IF (ASSOCIATED(p)) THEN
                    stdout("********************************")
                    stdout("patch     istart_p ",'p%istart_p(1:2)')
                    stdout("patch     iend_p ",'p%iend_p(1:2)')
                    stdout("********************************")
                ENDIF
            DO WHILE (ASSOCIATED(p))
                stdout("--------------------------------")
                stdout("patch     istart ",'p%istart(1:2)')
                stdout("patch     iend   ",'p%iend(1:2)')
                stdout("patch     g_start",'p%ghostsize(1)','p%ghostsize(3)')
                stdout("patch     g_end  ",'p%ghostsize(2)','p%ghostsize(4)')
                stdout("patch alloc from ",'p%lo_a(1:2)')
                stdout("            to   ",'p%hi_a(1:2)')
                stdout("--------------------------------")
                p => Mesh1%subpatch%next()
            ENDDO
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
        ENDIF

        !Fill in the allocated field arrays (incl. ghost nodes) with some data
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                Field1_n(1) = -10.0_MK * (rank+1)
                Field1_n(2) = -10.0_MK * (rank+1) - 1
                Field2_n    = -1.0_MK
        end foreach

        !Change the values of the real nodes to something that depends
        ! on position (so that we can later check the values)
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
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
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
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
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .LT. 0.0_MK) THEN
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
                Assert_Equal_Within(Field1_n(1) ,COS(2._MK*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,COS(2._MK*pi*pos(1)) + 2.0_MK,1e-5)
                Assert_Equal_Within(Field2_n    ,1.0_MK,1e-5)
        end foreach
        Assert_Equal(nb_errors,0)
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        CALL Mesh1%destroy(info)
        Assert_Equal(info,0)
        CALL Field1%destroy(info)
        Assert_Equal(info,0)
        CALL Field2%destroy(info)
        Assert_Equal(info,0)

        end_subroutine()
    end test


end test_suite
