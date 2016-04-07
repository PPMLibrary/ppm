test_suite ppm_module_neighlist
  USE ppm_module_interfaces

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
  REAL(MK),PARAMETER              :: pi = 3.1415926535897931_MK
  REAL(MK),PARAMETER              :: skin = 0._MK
  INTEGER,PARAMETER               :: ndim=2
  INTEGER,PARAMETER               :: pdim=2
  INTEGER                         :: decomp,assig,tolexp
  REAL(MK)                        :: tol,min_rcp,max_rcp
  INTEGER                         :: info,comm,rank,nproc
  INTEGER                         :: topoid
  INTEGER                         :: np = 100000
  INTEGER                         :: npart
  INTEGER                         :: mp
  INTEGER                         :: newnp
  REAL(MK),DIMENSION(:,:),POINTER :: xp => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: rcp => NULL()
  REAL(MK),DIMENSION(:,:),POINTER :: wp => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: h => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: p_h => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: len_phys => NULL()
  REAL(MK),DIMENSION(:  ),POINTER :: ghostlayer => NULL()
  INTEGER, DIMENSION(:  ),POINTER :: ghostsize => NULL()
  INTEGER                         :: i,j,k,sum1,sum2
  INTEGER                         :: p_i
  INTEGER, DIMENSION(6)           :: bcdef
  REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
  INTEGER, DIMENSION(:  ),POINTER :: nm => NULL()
  INTEGER                         :: seedsize
  INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
  REAL(MK), DIMENSION(:),ALLOCATABLE :: randnb
  INTEGER                          :: isymm = 0
  LOGICAL                          :: lsymm = .FALSE.,ok
  REAL(MK)                         :: t0,t1,t2,t3
  REAL(MK)                         :: eps

  init
    USE ppm_module_topo_typedef
    USE ppm_module_init

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
    &     ghostsize(ndim),ghostlayer(2*ndim),&
    &     nm(ndim),h(ndim),p_h(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    max_phys(1:ndim) = 1.0_MK
    len_phys(1:ndim) = max_phys-min_phys
    ghostsize(1:ndim) = 2
    ghostlayer(1:2*ndim) = max_rcp
    bcdef(1:6) = ppm_param_bcdef_periodic

    eps = EPSILON(1.0_MK)
    tolexp = INT(LOG10(EPSILON(1.0_MK)))

    NULLIFY(xp,rcp,wp)

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

  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)

    DEALLOCATE(min_phys,max_phys,len_phys,ghostsize,nm)
  end finalize


  setup
    CALL RANDOM_SEED(size=seedsize)
    ALLOCATE(seed(seedsize))
    ALLOCATE(randnb((1+ndim)*np),STAT=info)
    DO i=1,seedsize
        seed(i)=10+i*i*(rank+1)
    ENDDO
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(randnb)

    ALLOCATE(xp(ndim,np),rcp(np),wp(pdim,np),STAT=info)
  end setup


  teardown
    DEALLOCATE(xp,rcp,wp,STAT=info)
    DEALLOCATE(seed,randnb)
  end teardown

  test memleak
    USE ppm_module_topo_typedef
    USE ppm_module_MKtopo
    USE ppm_module_map
    USE ppm_module_topo_check
    USE ppm_module_util_dbg
    USE ppm_module_test

    INTEGER             :: mpart
    INTEGER             :: newnpart
    INTEGER             :: oldip,ip = -1
    REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
    REAL(MK),DIMENSION(ndim)    :: cp
    REAL(MK)            :: gl = 0.001_MK
    REAL(MK), PARAMETER         :: skin = 0.0_MK
    REAL(MK)            :: h
    INTEGER, DIMENSION(:),POINTER   :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER :: vlist => NULL()
    INTEGER             :: snpart
    INTEGER             :: i,j,k
    INTEGER, DIMENSION(:),POINTER   :: pidx => NULL()
    TYPE(ppm_t_clist),DIMENSION(:),POINTER :: clist => NULL()

    npart=1000
    CALL part_init(p,npart,min_phys,max_phys,info,&
    &    ppm_param_part_init_cartesian,0.5_MK)
    h = 2.0_MK*(len_phys(1)/(SQRT(REAL(npart,mk))))
    gl = 0.0_MK
    bcdef(1:6) = ppm_param_bcdef_freespace
    NULLIFY(nvlist,vlist,pidx)

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal

    topoid = 0

    CALL ppm_MKtopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef,gl+skin,cost,info)

    CALL ppm_map_part_global(topoid,p,npart,info)
    CALL ppm_map_part_send(npart,newnpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
    npart=newnpart


    CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl+skin,info)
    CALL ppm_map_part_send(npart,mpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

    CALL ppm_topo_check(topoid,p,npart,ok,info)
    !CALL ppm_dbg_print(topoid,gl+skin,1,1,info,p,npart,mpart)

    ALLOCATE(pidx(npart))
    FORALL(k=1:npart) pidx(k) = k

    CALL ppm_neighlist_vlist(topoid,p,mpart,h,skin,.TRUE.,vlist,nvlist,info,pidx,npart,clist)
    Assert_Equal(info,0)

    CALL ppm_clist_destroy(clist,info)
    Assert_Equal(info,0)

  end test

  !uncomment for longer, more thorough testing.
  !test stack_overflow({npart: [10,1000,10000,100000,1000000,10000000]})
  test stack_overflow({npart: [10,1000,10000,100000]})
    USE ppm_module_topo_typedef
    USE ppm_module_MKtopo
    USE ppm_module_map
    USE ppm_module_topo_check
    USE ppm_module_util_dbg
    USE ppm_module_test

    INTEGER             :: mpart
    INTEGER             :: newnpart
    INTEGER             :: oldip,ip = -1
    REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
    REAL(MK),DIMENSION(ndim)    :: cp
    REAL(MK)            :: gl = 0.001_MK
    REAL(MK), PARAMETER         :: skin = 0.0_MK
    REAL(MK)            :: h
    INTEGER, DIMENSION(:),POINTER   :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER :: vlist => NULL()
    INTEGER             :: snpart
    INTEGER             :: i,j,k
    INTEGER, DIMENSION(:),POINTER   :: pidx

    CALL part_init(p,npart,min_phys,max_phys,info,&
    &    ppm_param_part_init_cartesian,0.5_MK)
    !print *,npart
    h = 2.0_MK*(len_phys(1)/(SQRT(REAL(npart,mk))))
    gl = 0.0_MK
    bcdef(1:6) = ppm_param_bcdef_freespace
    NULLIFY(nvlist,vlist,pidx)

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal

    topoid = 0

    CALL ppm_MKtopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef,gl+skin,cost,info)

    CALL ppm_map_part_global(topoid,p,npart,info)
    CALL ppm_map_part_send(npart,newnpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
    npart=newnpart


    CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl+skin,info)
    CALL ppm_map_part_send(npart,mpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

    CALL ppm_topo_check(topoid,p,npart,ok,info)
    !CALL ppm_dbg_print(topoid,gl+skin,1,1,info,p,npart,mpart)

    ALLOCATE(pidx(npart))
    FORALL(k=1:npart) pidx(k) = k

    CALL ppm_neighlist_vlist(topoid,p,mpart,h,skin,.TRUE.,&
    &            vlist,nvlist,info)!,pidx)

    Assert_Equal(info,0)
    deALLOCATE(p,vlist,nvlist,pidx)
  end test

  test symbcvlistsize
    ! tests symmetric boundary conditions and vlist size/content
    USE ppm_module_topo_typedef
    USE ppm_module_MKtopo
    USE ppm_module_topo_check
    USE ppm_module_util_dbg
    USE ppm_module_map
    !  INTEGER             :: npart = 20**2
    INTEGER             :: npart = 20
    INTEGER             :: newnpart
    INTEGER             :: mpart
    REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
    REAL(MK),DIMENSION(  :),POINTER :: w => NULL()
    REAL(MK), PARAMETER         :: gl = 0.1_MK
    REAL(MK)            :: h
    REAL(MK), PARAMETER         :: skin = 0.01_MK
    INTEGER, DIMENSION(:),POINTER   :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER :: vlist => NULL()

    ALLOCATE(p(2,npart))
    h = 0.05_MK
    DO i=1,npart
        p(1,i) = h/2.0_MK + h*(i-1)
        p(2,i) = h/2.0_MK
    ENDDO

    bcdef(1:2) = ppm_param_bcdef_freespace
    bcdef(3:4) = ppm_param_bcdef_symmetry
    bcdef(5:6) = ppm_param_bcdef_freespace

    ALLOCATE(w(npart))
    w(:) = rank+1
#ifdef __MPI
    CALL MPI_BARRIER(comm,info)
#endif
    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_ypencil
    !decomp = ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal

    topoid = 0

    CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef,gl,cost,info)

    CALL ppm_map_part_global(topoid,p,npart,info)
    CALL ppm_map_part_push(w,npart,info)
    CALL ppm_map_part_send(npart,newnpart,info)
    CALL ppm_map_part_pop(w,npart,newnpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
    npart=newnpart

    CALL ppm_topo_check(topoid,p,npart,ok,info)

    Assert_True(ok)
    !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart)

    CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,1,gl,info)
    CALL ppm_map_part_push(w,npart,info)
    CALL ppm_map_part_send(npart,mpart,info)
    CALL ppm_map_part_pop(w,npart,mpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

    CALL ppm_topo_check(topoid,p,npart,ok,info)

    Assert_True(ok)
    CALL ppm_neighlist_vlist(topoid,p,mpart,gl/2.0_MK,skin,.TRUE.,&
    &            vlist,nvlist,info)

    !CALL ppm_dbg_print(topoid,gl,1,nvlist,info,p,npart,mpart)

    Assert_Equal(vlist(1,1),21)
    IF (nproc.EQ.2) THEN
       !Assert_Equal(vlist(2,1),31)
       !Assert_Equal(vlist(2,10),40)
       !Assert_Equal(vlist(1,20),30)
    ENDIF
    Assert_Equal(vlist(1,10),30)
  end test

  test symBC_neighlist
    ! tests symmetric boundary conditions and ghost get
    USE ppm_module_topo_typedef
    USE ppm_module_MKtopo
    USE ppm_module_map
    USE ppm_module_topo_check
    USE ppm_module_util_dbg

    INTEGER             :: npart = 4
    INTEGER             :: newnpart
    INTEGER             :: mpart
    INTEGER             :: oldip,ip = -1
    REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
    REAL(MK),DIMENSION(ndim)    :: cp
    REAL(MK), PARAMETER         :: gl = 0.1_MK
    REAL(MK), PARAMETER         :: skin = 0.05_MK
    INTEGER, DIMENSION(:),POINTER   :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER :: vlist => NULL()

    ALLOCATE(p(ndim,npart))
    p(1,1) = 0.05_MK
    p(2,1) = 0.5_MK
    p(1,2) = 0.5_MK
    p(2,2) = 0.05_MK
    p(1,3) = 0.05_MK
    p(2,3) = 0.05_MK
    p(1,4) = 0.5_MK
    p(2,4) = 0.95_MK

    bcdef(1:6) = ppm_param_bcdef_symmetry
    NULLIFY(nvlist,vlist)

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_cuboid
    !decomp = ppm_param_decomp_xpencil
    assig  = ppm_param_assign_internal

    topoid = 0

    CALL ppm_MKtopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
    &           gl+skin,cost,info)

    CALL ppm_map_part_global(topoid,p,npart,info)
    CALL ppm_map_part_send(npart,newnpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
    npart=newnpart


    CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl+skin,info)
    CALL ppm_map_part_send(npart,mpart,info)
    CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

    CALL ppm_topo_check(topoid,p,npart,ok,info)
    !CALL ppm_dbg_print(topoid,gl+skin,1,1,info,p,npart,mpart)

    CALL ppm_neighlist_vlist(topoid,p,mpart,gl,skin,.TRUE.,&
    &            vlist,nvlist,info)

    !DO i=1,mpart
    !    print *,i,p(:,i)
    !    DO j=1,nvlist(i)
    !    print *,'    ',vlist(j,i)
    !    ENDDO
    !ENDDO

    ! DO the tests
    ! p(1)
    cp(1) = -0.05_MK
    cp(2) = 0.5_MK
    ip = -1
    DO i=npart+1,mpart
        IF ((ABS(p(1,i)-cp(1)).lt.eps).and.&
        &   (ABS(p(2,i)-cp(2)).lt.eps)) THEN
          ip = i
          exit
        ENDIF
    ENDDO
    IF (nproc.eq.1) THEN
       Assert_False(ip.eq.-1)
       Assert_True((vlist(1,1).eq.ip).or.(vlist(1,ip).eq.1))
    ENDIF
    ! p(2)
    cp(1) = 0.5_MK
    cp(2) = -0.05_MK
    ip = -1
    DO i=npart+1,mpart
        IF ((ABS(p(1,i)-cp(1)).lt.eps).and.&
        &   (ABS(p(2,i)-cp(2)).lt.eps)) THEN
        ip = i
        exit
        ENDIF
    ENDDO
    IF (nproc.eq.1) THEN
       Assert_False(ip.eq.-1)
       Assert_True((vlist(1,2).eq.ip).or.(vlist(1,ip).eq.2))
    ENDIF
    ! p(3)
    cp(1) = 0.05_MK
    cp(2) = -0.05_MK
    ip = -1
    DO i=npart+1,mpart
        IF ((ABS(p(1,i)-cp(1)).lt.eps).and.&
        &   (ABS(p(2,i)-cp(2)).lt.eps)) THEN
        ip = i
        exit
        ENDIF
    ENDDO
    IF (nproc.eq.1) THEN
       Assert_False(ip.eq.-1)
       ok = (vlist(1,3).eq.ip).or.&
       &    (vlist(2,3).eq.ip).or.&
       &    (vlist(1,ip).eq.3).or.&
       &    (vlist(2,ip).eq.3)
       Assert_True(ok)
    ENDIF

    cp(1) = -0.05_MK
    cp(2) = 0.05_MK
    ip = -1
    DO i=npart+1,mpart
        IF ((ABS(p(1,i)-cp(1)).lt.eps).and.&
        &   (ABS(p(2,i)-cp(2)).lt.eps)) THEN
        ip = i
        exit
        ENDIF
    ENDDO
    IF (nproc.eq.1) THEN
       Assert_False(ip.eq.-1)
       ok = (vlist(1,3).eq.ip).or.&
       &    (vlist(2,3).eq.ip).or.&
       &    (vlist(1,ip).eq.3).or.&
       &    (vlist(2,ip).eq.3)
       Assert_True(ok)
    ENDIF
    cp(1) = -0.05_MK
    cp(2) = -0.05_MK
    ip = -1
    DO i=npart+1,mpart
        IF ((ABS(p(1,i)-cp(1)).lt.eps).and.&
        &   (ABS(p(2,i)-cp(2)).lt.eps)) THEN
        ip = i
        exit
        ENDIF
    ENDDO
    IF (nproc.eq.1) THEN
       Assert_False(ip.eq.-1)
       ok = (vlist(1,3).eq.ip).or.&
       &    (vlist(2,3).eq.ip).or.&
       &    (vlist(1,ip).eq.3).or.&
       &    (vlist(2,ip).eq.3)
       Assert_True(ok)
    ENDIF

    ! p(4)
    cp(1) = 0.5_MK
    cp(2) = 1.05_MK
    ip = -1
    DO i=npart+1,mpart
        IF ((ABS(p(1,i)-cp(1)).lt.eps).and.&
        &   (ABS(p(2,i)-cp(2)).lt.eps)) THEN
        ip = i
        exit
        ENDIF
    ENDDO
    IF (nproc.eq.1) THEN
       Assert_False(ip.eq.-1)
       Assert_True((vlist(1,4).eq.ip).or.(vlist(1,ip).eq.4))
    ENDIF
  end test

end test_suite