test_suite ppm_module_particles

  USE ppm_module_particles_typedef
  USE ppm_module_topo_typedef
  USE ppm_module_field_typedef
  USE ppm_module_operator_typedef
  USE ppm_module_interfaces
  USE ppm_module_data

  INTEGER, PARAMETER :: MK = KIND(1.0d0) !KIND(1.0e0)

  COMPLEX(MK), DIMENSION(:),   POINTER :: wp_1c => NULL()
  COMPLEX(MK), DIMENSION(:,:), POINTER :: wp_2c => NULL()

  REAL(MK),              PARAMETER                   :: tol=EPSILON(1._MK)*100
  REAL(MK),              PARAMETER                   :: pi = ACOS(-1._MK)
  REAL(MK),              PARAMETER                   :: skin = 0._MK
  REAL(MK)                                           :: cutoff
  REAL(MK),              DIMENSION(:,:), POINTER     :: xp=>NULL()
  REAL(MK),              DIMENSION(:  ), POINTER     :: min_phys=>NULL(),max_phys=>NULL()
  REAL(MK),              DIMENSION(:  ), POINTER     :: len_phys=>NULL()
  REAL(MK),              DIMENSION(:  ), POINTER     :: cost=>NULL()
  REAL(MK)                                           :: t0,t1,t2,t3
  REAL(MK)                                           :: err
  REAL(MK),              DIMENSION(:),   POINTER     :: wp_1r => NULL()
  REAL(MK),              DIMENSION(:,:), POINTER     :: wp_2r => NULL()
  REAL(ppm_kind_double), DIMENSION(:),   ALLOCATABLE :: coeffs

  INTEGER,                 PARAMETER                   :: debug = 0
  INTEGER,                 PARAMETER                   :: ndim=2
  INTEGER                                              :: decomp,assig,tolexp
  INTEGER                                              :: info,comm,rank,nproc,topoid
  INTEGER                                              :: np_global = 3000
  INTEGER                                              :: i,j,k,ip,wp_id
  INTEGER                                              :: nstep
  INTEGER,                 DIMENSION(3)                :: ldc
  INTEGER,                 DIMENSION(6)                :: bcdef
  INTEGER                                              :: isymm = 0
  INTEGER                                              :: seedsize
  INTEGER,                 DIMENSION(:),   ALLOCATABLE :: seed
  INTEGER,                 DIMENSION(:),   POINTER     :: nvlist=>NULL()
  INTEGER,                 DIMENSION(:,:), POINTER     :: vlist=>NULL()
  INTEGER,                 DIMENSION(:),   POINTER     :: wp_1i => NULL()
  INTEGER,                 DIMENSION(:,:), POINTER     :: wp_2i => NULL()
  INTEGER(ppm_kind_int64), DIMENSION(:),   POINTER     :: wp_1li => NULL()
  INTEGER(ppm_kind_int64), DIMENSION(:,:), POINTER     :: wp_2li => NULL()
  INTEGER,                 DIMENSION(:),   ALLOCATABLE :: degree,order
  INTEGER                                              :: nterms

  LOGICAL, DIMENSION(:), POINTER :: wp_1l => NULL()

  init
    USE ppm_module_init
    USE ppm_module_mktopo

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    max_phys(1:ndim) = 1.0_MK
    max_phys(ndim) = 1.4_MK
    len_phys(1:ndim) = max_phys-min_phys
    bcdef(1:6) = ppm_param_bcdef_periodic

#ifdef __MPI
    comm = MPI_COMM_WORLD
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    rank = 0
    nproc = 1
#endif
    tolexp = INT(LOG10(EPSILON(1._MK)))+10
    CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

    CALL RANDOM_SEED(SIZE=seedsize)
    ALLOCATE(seed(seedsize))
    DO i=1,seedsize
       seed(i)=10+i*i*(rank+1)
    ENDDO
    CALL RANDOM_SEED(put=seed)

    !----------------
    ! make topology
    !----------------
    decomp =  ppm_param_decomp_cuboid
    !decomp = ppm_param_decomp_xpencil
    assig  = ppm_param_assign_internal

    topoid = 0
    IF (ndim.EQ.2) THEN
       cutoff = 0.15_MK
    ELSE
       cutoff = 0.25_MK
    ENDIF

    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
  end init


  finalize
    USE ppm_module_finalize
    CALL ppm_finalize(info)
    DEALLOCATE(min_phys,max_phys,len_phys)
  end finalize

  setup
  end setup

  teardown
  end teardown

  test ghost_mappings
    USE ppm_module_io_vtk

    TYPE(ppm_t_particles_d), TARGET :: Part1

    TYPE(ppm_t_field) :: Field1
    TYPE(ppm_t_field) :: Field2
    TYPE(ppm_t_field) :: Field3

    CLASS(ppm_t_neighlist_d_), POINTER :: Nlist => NULL()

    CLASS(ppm_t_discr_data), POINTER :: prop => NULL()

    start_subroutine("ghost_mappings")

    !--------------------------
    !Define Fields
    !--------------------------
    CALL Field1%create(5,info,name="VecF1") !vector field
    Assert_Equal(info,0)
    CALL Field2%create(1,info,name="ScaF1") !scalar field
    Assert_Equal(info,0)
    CALL Field3%create(3,info,name="VecF2") !vector field
    Assert_Equal(info,0)

    CALL Part1%initialize(np_global,info,topoid=topoid,name="Part1")
    Assert_Equal(info,0)

    CALL Part1%create_neighlist(Part1,info)
    Assert_Equal(info,0)

    CALL Part1%set_cutoff(3._MK * Part1%h_avg,info)
    Assert_Equal(info,0)

    ALLOCATE(wp_2r(ndim,Part1%Npart))
    CALL RANDOM_NUMBER(wp_2r)
    wp_2r = (wp_2r - 0.5_MK) * Part1%h_avg * 0.15_MK

    CALL Part1%move(wp_2r,info)
    Assert_Equal(info,0)
    DEALLOCATE(wp_2r)

    Assert_True(Part1%has_neighlist(Part1))
    CALL Part1%apply_bc(info)
    Assert_Equal(info,0)

    CALL Part1%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)

    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    CALL Field1%discretize_on(Part1,info)
    Assert_Equal(info,0)
    CALL Field2%discretize_on(Part1,info)
    Assert_Equal(info,0)
    CALL Field3%discretize_on(Part1,info)
    Assert_Equal(info,0)

    Assert_True(Part1%has_ghosts(Field1))
    Assert_True(Part1%has_ghosts(Field2))
    Assert_True(Part1%has_ghosts(Field3))

    CALL Part1%comp_neighlist(info)
    Assert_Equal(info,0)

    Nlist => Part1%get_neighlist(Part1)
    Assert_True(ASSOCIATED(Nlist))
    Nlist => NULL()

    !Set field values
    foreach p in particles(Part1) with positions(x) sca_fields(F2=Field2) vec_fields(F1=Field1,F3=Field3)
        F1_p(1) = x_p(1)
        F1_p(2) = x_p(2)
        F1_p(3) = -3._MK
        F1_p(4) = -4._MK
        F1_p(5) = -4._MK
        F2_p    = -10._MK
        F3_p(1) = -11._MK
        F3_p(2) = -12._MK
        F3_p(3) = -13._MK
    end foreach

    !  print particles to a VTK file
    !CALL ppm_vtk_particles("output",Part1,info)
    !Assert_Equal(info,0)


    !Check that PPM knows that the ghosts are now invalid for all the fields
    Assert_False(Part1%has_ghosts(Field1))
    Assert_False(Part1%has_ghosts(Field2))
    Assert_False(Part1%has_ghosts(Field3))


    ! Do a ghost mapping, but only for fields 2 and 3.
    CALL Part1%map_ghost_get(info)
    Assert_Equal(info,0)
    CALL Part1%map_ghost_push(info,Field2)
    Assert_Equal(info,0)
    CALL Part1%map_ghost_push(info,Field3)
    Assert_Equal(info,0)

    CALL Part1%map_ghost_send(info)
    Assert_Equal(info,0)

    CALL Part1%map_ghost_pop(info,Field3)
    Assert_Equal(info,0)
    CALL Part1%map_ghost_pop(info,Field2)
    Assert_Equal(info,0)

    CALL Part1%map_ghost_pop_positions(info)
    Assert_Equal(info,0)

    ! Check the states (ghosts should be ok for Field2 and Field3
    ! but not for Field1)
    Assert_False(Part1%has_ghosts(Field1))
    Assert_True(Part1%has_ghosts(Field2))
    Assert_True(Part1%has_ghosts(Field3))

    !Move particles
    ALLOCATE(wp_2r(ndim,Part1%Npart))
    CALL RANDOM_NUMBER(wp_2r)
    wp_2r = (wp_2r - 0.5_MK) * Part1%ghostlayer * 0.99_MK
    CALL Part1%move(wp_2r,info)
    Assert_Equal(info,0)
    DEALLOCATE(wp_2r)
    CALL Part1%apply_bc(info)
    Assert_Equal(info,0)

    ! Do a partial mapping, but only map fields 2 and 3.
    CALL Part1%map_positions(info)
    Assert_Equal(info,0)

    CALL Part1%map_push(info,Field2)
    Assert_Equal(info,0)
    CALL Part1%map_push(info,Field3)
    Assert_Equal(info,0)

    CALL Part1%map_send(info)
    Assert_Equal(info,0)

    !NON-Blocking
    !CALL Part1%map_isend(info)
    !Assert_Equal(info,0)


    CALL Part1%map_pop(info,Field3)
    Assert_Equal(info,0)
    CALL Part1%map_pop(info,Field2)
    Assert_Equal(info,0)

    CALL Part1%map_pop_positions(info)
    Assert_Equal(info,0)

    !Check that are still correct
    foreach p in particles(Part1) with positions(x) sca_fields(F2=Field2) vec_fields(F3=Field3)
        Assert_Equal_Within(F2_p   , -10._MK,tol)
        Assert_Equal_Within(F3_p(1), -11._MK,tol)
        Assert_Equal_Within(F3_p(2), -12._MK,tol)
        Assert_Equal_Within(F3_p(3), -13._MK,tol)
    end foreach

    !CALL Part1%destroy(info)
    !Assert_Equal(info,0)
    CALL Field1%destroy(info)
    Assert_Equal(info,0)
    CALL Field2%destroy(info)
    Assert_Equal(info,0)
    CALL Field3%destroy(info)
    Assert_Equal(info,0)

    CALL Part1%destroy(info)
    Assert_Equal(info,0)

    end_subroutine()
  end test

  test part_growth_and_shrink
    TYPE(ppm_t_particles_d), TARGET :: Part1

    TYPE(ppm_t_field) :: Field1
    TYPE(ppm_t_field) :: Field2

    CLASS(ppm_t_discr_data), POINTER :: prop => NULL()

    INTEGER, DIMENSION(:), POINTER :: list_del_parts

    REAL(MK) :: tmp

    start_subroutine("part_growth_and_shrink")

    !--------------------------
    ! Define vector Fields
    !--------------------------
    CALL Field1%create(5,info,name="VecF1") !vector field
    Assert_Equal(info,0)
    CALL Field2%create(3,info,name="VecF2") !vector field
    Assert_Equal(info,0)

    np_global=np_global*2

    CALL Part1%initialize(np_global,info,topoid=topoid,name="Part1")
    Assert_Equal(info,0)

    CALL Part1%apply_bc(info)
    Assert_Equal(info,0)

    Part1%ghostlayer=2._MK*Part1%h_avg

    CALL Part1%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)


    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    CALL Field1%discretize_on(Part1,info)
    Assert_Equal(info,0)
    CALL Field2%discretize_on(Part1,info)
    Assert_Equal(info,0)

    Assert_True(Part1%has_ghosts(Field1))
    Assert_True(Part1%has_ghosts(Field2))

    !Set field values
    foreach p in particles(Part1) with positions(x) vec_fields(F1=Field1,F2=Field2)
        F1_p(1) = x_p(1)
        F1_p(2) = x_p(2)
        F1_p(3) = -3._MK
        F1_p(4) = -4._MK
        F1_p(5) = -4._MK
        F2_p(1) = -11._MK
        F2_p(2) = -12._MK
        F2_p(3) = -13._MK
    end foreach

    !Check that PPM knows that the ghosts are now invalid for all the fields
    Assert_False(Part1%has_ghosts(Field1))
    Assert_False(Part1%has_ghosts(Field2))

    ! Do a ghost mapping, but only for fields 2.
    CALL Part1%map_ghost_get(info)
    Assert_Equal(info,0)
    CALL Part1%map_ghost_push(info,Field2)
    Assert_Equal(info,0)

    CALL Part1%map_ghost_send(info)
    Assert_Equal(info,0)

    CALL Part1%map_ghost_pop(info,Field2)
    Assert_Equal(info,0)
    CALL Part1%map_ghost_pop_positions(info)
    Assert_Equal(info,0)

    ! Check the states (ghosts should be ok for Field2 and Field3
    ! but not for Field1)
    Assert_False(Part1%has_ghosts(Field1))
    Assert_True(Part1%has_ghosts(Field2))

    !Check that are still correct
    foreach p in particles(Part1) with positions(x) vec_fields(F2=Field2)
        Assert_Equal_Within(F2_p(1), -11._MK,tol)
        Assert_Equal_Within(F2_p(2), -12._MK,tol)
        Assert_Equal_Within(F2_p(3), -13._MK,tol)
    end foreach

    stdout("Npart=",'Part1%Npart'," Particle array size=",'Part1%size()')

    CALL Part1%grow_size(info)
    Assert_Equal(info,0)

    stdout("Npart after growth=",'Part1%Npart'," Particle array size=",'Part1%size()')

    CALL Part1%del_parts(Part1%Npart,info)
    Assert_Equal(info,0)

    stdout("Npart after 1 particle delete, is=",'Part1%Npart'," Particle array size=",'Part1%size()')

    ALLOCATE(list_del_parts(Part1%Npart),STAT=info)
    Assert_Equal(info,0)

    j=0
    DO i=1,Part1%Npart
       CALL RANDOM_NUMBER(tmp)
       IF (tmp.GT.0.6_MK) THEN
          j=j+1
          k=MOD(INT(tmp*10000),Part1%Npart)+1
          list_del_parts(j)=k
       ENDIF
    ENDDO

    stdout("Number of particles (with repetition) to be deleted are =",j)

    CALL Part1%del_parts(list_del_parts,j,info)
    Assert_Equal(info,0)

    stdout("Npart after",j," particles delete, is=",'Part1%Npart'," Particle array size=",'Part1%size()')

    DEALLOCATE(list_del_parts,STAT=info)
    Assert_Equal(info,0)

    !--------------------------
    ! Destro Field
    !--------------------------
    CALL Field1%destroy(info)
    Assert_Equal(info,0)
    CALL Field2%destroy(info)
    Assert_Equal(info,0)
    !--------------------------
    ! Destro particle
    !--------------------------
    CALL Part1%destroy(info)
    Assert_Equal(info,0)

    end_subroutine()
  end test

  !-------------------------------------------------------------
  ! test function
  !-------------------------------------------------------------
  PURE FUNCTION f0_test(pos,ndim)
  IMPLICIT NONE
  INTEGER,                   INTENT(IN   ) :: ndim
  REAL(MK), DIMENSION(ndim), INTENT(IN   ) :: pos
  REAL(MK)                                 :: f0_test
  f0_test=SIN(2._MK*pi*pos(1))*COS(2._MK*pi*pos(2))*SIN(2._MK*pi*pos(ndim))
  RETURN
  END FUNCTION f0_test

end test_suite
