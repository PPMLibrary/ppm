test_suite ppm_module_mesh_small

USE ppm_module_mesh_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_particles_typedef
USE ppm_module_mktopo

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: mk = kind(1.0d0) !kind(1.0e0)
REAL(MK),PARAMETER              :: pi = ACOS(-1._mk)
INTEGER,PARAMETER               :: ndim=2
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc
REAL(MK)                        :: tol
INTEGER                         :: topoid=-1
REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()

INTEGER, DIMENSION(:  ),POINTER :: ighostsize => NULL()
REAL(MK)                        :: sca_ghostsize
INTEGER                         :: seedsize
INTEGER,  DIMENSION(:),allocatable :: seed

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
CLASS(ppm_t_subpatch_),POINTER   :: patch => NULL()

INTEGER                          :: mypatchid
REAL(MK),DIMENSION(6)       :: my_patch
REAL(MK),DIMENSION(ndim)         :: offset

INTEGER, DIMENSION(:), POINTER                 :: wp_1i => NULL()
INTEGER, DIMENSION(:,:), POINTER               :: wp_2i => NULL()
INTEGER(ppm_kind_int64),DIMENSION(:),  POINTER :: wp_1li => NULL()
INTEGER(ppm_kind_int64),DIMENSION(:,:),POINTER :: wp_2li => NULL()
REAL(MK), DIMENSION(:),   POINTER              :: wp_1r => NULL()
REAL(MK), DIMENSION(:,:), POINTER              :: wp_2r => NULL()
complex(mk), DIMENSION(:),   POINTER           :: wp_1c => NULL()
complex(mk), DIMENSION(:,:), POINTER           :: wp_2c => NULL()
LOGICAL, DIMENSION(:),   POINTER               :: wp_1l => NULL()

!---------------- init -----------------------

    init

        USE ppm_module_topo_typedef
        USE ppm_module_init

        ALLOCATE(min_phys(ndim),max_phys(ndim),&
            &         ighostsize(ndim),nm(ndim),h(ndim))

        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
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

        CALL RANDOM_SEED(size=seedsize)
        ALLOCATE(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        CALL RANDOM_SEED(put=seed)

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

    test small_test
        TYPE(ppm_t_field) :: Field1,Field2
        TYPE(ppm_t_particles_d),TARGET :: Part1

        REAL(ppm_kind_double),DIMENSION(ndim) :: x0

        REAL(ppm_kind_double),DIMENSION(ndim) :: pos
        REAL(ppm_kind_double),DIMENSION(ndim) :: cutoff
        INTEGER :: p_idx, nb_errors
        INTEGER :: np_global = 3000000
        CLASS(ppm_t_discr_info_),POINTER              :: dinfo => NULL()
        REAL(ppm_kind_double),DIMENSION(:  ), POINTER :: up_1d => NULL()
        REAL(ppm_kind_double),DIMENSION(:,:), POINTER :: up_2d => NULL()
        LOGICAL :: assoc

        start_subroutine("test_interp")

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.07_mk
        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
        &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)

        Nm = 35
        Nm(ndim) = 65
        CALL Mesh1%create(topoid,offset,info,Nm=Nm,&
        &   ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        CALL Part1%initialize(np_global,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

        CALL Part1%create_neighlist(Part1,info)
        or_fail("failed to create neighbour list")

        CALL Part1%set_cutoff(3._mk * Part1%h_avg,info)
        Assert_Equal(info,0)

        ALLOCATE(wp_2r(ndim,Part1%Npart))
        CALL RANDOM_NUMBER(wp_2r)
        wp_2r = (wp_2r - 0.5_mk) * Part1%h_avg * 0.15_mk
        CALL Part1%move(wp_2r,info)
        Assert_Equal(info,0)
        DEALLOCATE(wp_2r)

        Assert_true(Part1%has_neighlist(Part1))
        CALL Part1%apply_bc(info)
        Assert_Equal(info,0)

        CALL Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        CALL Part1%map_ghosts(info)
        Assert_Equal(info,0)


        if (ndim.eq.2) then
            my_patch(1:2*ndim) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.51_mk,0.99_mk,0.7_mk,0.78_mk/)
        endif
        CALL Mesh1%def_patch(my_patch,info)
        Assert_Equal(info,0)

        CALL Field1%create(3,info,name='vecField')
        Assert_Equal(info,0)
        CALL Field1%discretize_on(Mesh1,info)
        Assert_Equal(info,0)
        CALL Field2%create(1,info,name='scaField')
        Assert_Equal(info,0)
        CALL Field2%discretize_on(Mesh1,info)
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
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
            topo => ppm_topo(topoid)%t

            stdout("NB subdomains =  ",topo%nsubs)
            do i = 1,topo%nsublist
               isub = topo%isublist(i)
               stdout("coordinates subs Min =  ",'topo%min_subd(1:ndim,isub)')
               stdout("coordinates subs Max =  ",'topo%max_subd(1:ndim,isub)')
            enddo
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
            stdout("NB patch =  ",Mesh1%npatch)
            stdout("NB subpatch =  ",Mesh1%subpatch%nb)
            p => Mesh1%subpatch%begin()
            if (associated(p)) then
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
#ifdef __MPI
            CALL MPI_BARRIER(comm,info)
#endif
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
            !for real_and_ghosts
                pos(1:ndim) = sbpitr%get_pos(i,j)
                IF (Field2_n .lt. 0._mk) then
                    nb_errors = nb_errors + 1
                ENDIF
                Assert_Equal_Within(Field1_n(1) ,cos(2._mk*pi*pos(1)),        1e-5)
                Assert_Equal_Within(Field1_n(2) ,cos(2._mk*pi*pos(1)) + 2._mk,1e-5)
                Assert_Equal_Within(Field2_n    ,1._mk,1e-5)
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
        CALL Part1%destroy(info)
        Assert_Equal(info,0)

        end_subroutine()
    end test

!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function test_constant(pos,ndim) RESULT(res)
    REAL(MK)                              :: res
    INTEGER                 ,  intent(in) :: ndim
    REAL(MK), DIMENSION(ndim), intent(in) :: pos

    res =  1._mk !42.17_mk
end function

pure function test_linear(pos,ndim) RESULT(res)
    REAL(MK)                              :: res
    INTEGER                 ,  intent(in) :: ndim
    REAL(MK), DIMENSION(ndim), intent(in) :: pos

    res =  pos(1) + 10._mk*pos(2) + 100._mk*pos(ndim)
end function

pure function test_quadratic(pos,ndim) RESULT(res)
    REAL(MK)                              :: res
    INTEGER                 ,  intent(in) :: ndim
    REAL(MK), DIMENSION(ndim), intent(in) :: pos

    res =  pos(1)**2 + 10._mk*pos(2)**2 + 100._mk*pos(ndim)**2
end function

!!! check whether a particle is within a patch and more than a cutoff
!!! distance away from its boundaries.
pure function is_well_within(pos,patch,cutoff,ndim) RESULT(res)
    LOGICAL                               :: res
    REAL(MK), DIMENSION(ndim), intent(in) :: pos
    REAL(MK), DIMENSION(2*ndim),intent(in):: patch
    REAL(MK), DIMENSION(ndim), intent(in) :: cutoff
    INTEGER                 ,  intent(in) :: ndim

    res = ALL(pos(1:ndim).GE.(patch(1:ndim)+cutoff(1:ndim)))
    res = res .AND. ALL(pos(1:ndim).LE.(patch(ndim+1:2*ndim)-cutoff(1:ndim)))

end function

end test_suite
