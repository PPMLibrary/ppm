test_suite ppm_module_map_part

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: mk = kind(1.0d0) !kind(1.0e0)
  REAL(MK),PARAMETER              :: pi = 3.1415926535897931_mk
  REAL(MK),PARAMETER              :: skin = 0._mk
  INTEGER,PARAMETER               :: ndim=2
  INTEGER,PARAMETER               :: pdim=2
  INTEGER                         :: decomp,assig,tolexp
  REAL(MK)                        :: tol
  REAL(MK)                        :: min_rcp = 0.01_mk
  REAL(MK)                        :: max_rcp = 0.1_mk
  INTEGER                         :: info,comm,rank,nproc
  INTEGER                         :: topoid
  INTEGER, PARAMETER              :: np_init = 10000
  INTEGER                         :: np,mp
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
  INTEGER,  DIMENSION(:),allocatable :: seed
  REAL(MK), DIMENSION(:),allocatable :: randnb
  INTEGER                          :: isymm = 0
  LOGICAL                          :: lsymm = .false.,ok
  REAL(MK)                         :: t0,t1,t2,t3

    init

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_init

        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
        &         ghostsize(ndim),ghostlayer(2*ndim),&
        &         nm(ndim),h(ndim),p_h(ndim),STAT=info)

        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = max_rcp
        bcdef(1:6) = ppm_param_bcdef_periodic
        tol = EPSILON(1.0_mk)
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

        np = np_init
        CALL RANDOM_SEED(size=seedsize)
        ALLOCATE(seed(seedsize))
        ALLOCATE(randnb((1+ndim)*np),STAT=info)
        do i=1,seedsize
           seed(i)=10+i*i*(rank+1)
        enddo
        CALL RANDOM_SEED(put=seed)
        CALL RANDOM_NUMBER(randnb)

        ALLOCATE(xp(ndim,np),rcp(np),wp(pdim,np),STAT=info)

    end setup


    teardown

        DEALLOCATE(xp,rcp,wp,STAT=info)
        DEALLOCATE(seed,randnb)
        IF (associated(cost)) DEALLOCATE(cost)
        cost => NULL()

    end teardown

    test global
        ! test global mapping

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check

        !----------------
        ! create particles
        !----------------

        xp = 0.0_mk
        rcp = 0.0_mk

        !p_h = len_phys / real(npgrid,mk)
        !do j=1,npgrid
        !    do i=1,npgrid
        !        p_i = i + (j-1)*npgrid
        !        xp(1,p_i) = min_phys(1)+real(i-1,mk)*p_h(1)
        !        xp(2,p_i) = min_phys(2)+real(j-1,mk)*p_h(2)
        !        rcp(p_i) = min_rcp + (max_rcp-min_rcp)*randnb(p_i)
        !        do k=1,pdim
        !            wp(k,i) = rcp(i)*REAL(k,MK)
        !        enddo
        !    enddo
        !enddo
        do i=1,np
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo

        CALL RANDOM_NUMBER(xp)
        CALL RANDOM_NUMBER(wp)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)

        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(rcp,np,info)
        CALL ppm_map_part_push(wp,pdim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(wp,pdim,np,newnp,info)
        CALL ppm_map_part_pop(rcp,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        CALL ppm_topo_check(topoid,xp,newnp,ok,info)

        assert_true(ok)

    end test

    test partial
        ! test partial mapping

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check

        !----------------
        ! create particles
        !----------------

        xp = 0.0_mk
        rcp = 0.0_mk

        do i=1,np
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0


        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)

        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(rcp,np,info)
        CALL ppm_map_part_push(wp,pdim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(wp,pdim,np,newnp,info)
        CALL ppm_map_part_pop(rcp,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)
        np=newnp

        ! move all particles
        DEALLOCATE(randnb)
        ALLOCATE(randnb((1+ndim)*np))
        CALL RANDOM_NUMBER(randnb)

        do i=1,np
            do j=1,ndim
                xp(j,i) = xp(j,i) + (len_phys(j)/REAL(nproc,mk))*&
                &         randnb((ndim+1)*i-(ndim-j))
            enddo
        enddo

        ! do local mapping
        CALL ppm_map_part_partial(topoid,xp,np,info)
        CALL ppm_map_part_push(rcp,np,info)
        CALL ppm_map_part_push(wp,pdim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(wp,pdim,np,newnp,info)
        CALL ppm_map_part_pop(rcp,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)
        np=newnp

        CALL ppm_topo_check(topoid,xp,np,ok,info)
        assert_true(ok)

    end test

    test ghost_get
        ! test ghost get and how it treats pbc
        ! Ticket #32
        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_util_dbg
          INTEGER                         :: npart = 1
          INTEGER                         :: newnpart
          INTEGER                         :: mpart
          REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
          REAL(MK), PARAMETER             :: gl = 0.1_mk

        if (nproc.GT.1) return

        ALLOCATE(p(ndim,npart))
        p(1,1) = 0.05_mk
        p(2,1) = 0.05_mk
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        CALL ppm_map_part_global(topoid,p,npart,info)
        CALL ppm_map_part_send(npart,newnpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart

        CALL ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart)

        CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,2,1,info,p,npart,mpart)

    end test

    test pbc_and_map_load
        ! test ghost get map load and periodic BC

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_util_dbg
          INTEGER                         :: npart = 8
          INTEGER                         :: newnpart
          INTEGER                         :: mpart
          REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
          REAL(MK),DIMENSION(:)  ,POINTER :: w => NULL()
          REAL(MK),DIMENSION(2)           :: check
          REAL(MK), PARAMETER             :: gl = 0.1_mk

        if (nproc.GT.1) return

        ALLOCATE(p(ndim,npart),w(npart))
        p(1,1) = 0.05_mk  ! left
        p(2,1) = 0.5_mk
        p(1,2) = 0.95_mk  ! right
        p(2,2) = 0.5_mk
        p(1,3) = 0.5_mk   ! bottom
        p(2,3) = 0.05_mk
        p(1,4) = 0.5_mk   ! top
        p(2,4) = 0.95_mk

        p(1,5) = 0.05_mk  ! left-bottom
        p(2,5) = 0.05_mk
        p(1,6) = 0.05_mk  ! left-top
        p(2,6) = 0.95_mk
        p(1,7) = 0.95_mk  ! right-bottom
        p(2,7) = 0.05_mk
        p(1,8) = 0.95_mk  ! right-top
        p(2,8) = 0.95_mk

        w(:) = 1.0_mk

        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        CALL ppm_map_part_global(topoid,p,npart,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,newnpart,info)
        CALL ppm_map_part_pop(w,npart,newnpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart
        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart)

        CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(w,npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,2,1,info,p,npart,mpart)
        CALL ppm_map_part_store(info)

        CALL ppm_map_part_ghost_put(topoid,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_ghost_pop(w,1,npart,mpart,info)

        CALL ppm_map_part_load(info)
        CALL ppm_map_part_push(p,ndim,npart,info,.TRUE.)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,3,1,info,p,npart,mpart)

        assert_equal(mpart-npart,4+3*4) ! check number of ghosts

        ! now go through all particles and try to find their ghosts

        check(1) = 1.05_mk  ! left
        check(2) = 0.5_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = -0.05_mk  ! right
        check(2) =  0.5_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = 0.5_mk   ! bottom
        check(2) = 1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) =  0.5_mk   ! top
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))


        check(1) = 1.05_mk  ! left-bottom
        check(2) = 0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = 0.05_mk
        check(2) = 1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = 1.05_mk
        check(2) = 1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) =  1.05_mk  ! left-top
        check(2) =  0.95_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  0.05_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  1.05_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = -0.05_mk  ! right-bottom
        check(2) =  0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  0.95_mk
        check(2) =  1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = -0.05_mk
        check(2) =  1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = -0.05_mk  ! right-top
        check(2) =  0.95_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  0.95_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = -0.05_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        DEALLOCATE(w)

    end test

    function found_ghost(ghosts,n,cp)
      REAL(MK),DIMENSION(2,n),intent(in) :: ghosts
      INTEGER                ,intent(in) :: n
      REAL(MK),DIMENSION(2)  ,intent(in) :: cp
    LOGICAL                 :: found_ghost

      INTEGER                 :: i
      INTEGER                 :: found

    found_ghost = .false.

    found = 0
    do i=1,n
        if((abs(ghosts(1,i) - cp(1)).lt.tol).and. &
        &  (abs(ghosts(2,i) - cp(2)).lt.tol)) then
            found = found + 1
        endif
    enddo
    if (found.eq.1) then
        found_ghost = .true.
    endif

    end function found_ghost

    test getsymbc
        ! tests symmetric boundary conditions and ghost get
        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_util_dbg
          INTEGER                         :: npart = 1
          INTEGER                         :: newnpart
          INTEGER                         :: mpart
          REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
          REAL(MK), PARAMETER             :: gl = 0.1_mk

        if (nproc.GT.1) return

        ALLOCATE(p(ndim,npart))
        p(1,1) = 0.05_mk
        p(2,1) = 0.05_mk
        bcdef(1:6) = ppm_param_bcdef_symmetry

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        CALL ppm_map_part_global(topoid,p,npart,info)
        CALL ppm_map_part_send(npart,newnpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart

        CALL ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart)

        CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,2,1,info,p,npart,mpart)

    end test

    test get_mixed_bc
        ! tests symmetric and periodic BC conditions and ghost get
        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_util_dbg
          INTEGER                         :: npart = 256
          INTEGER                         :: snpart = 16
          INTEGER                         :: newnpart
          INTEGER                         :: mpart
          REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
          REAL(MK), PARAMETER             :: gl = 0.1_mk
          REAL(MK)                        :: h

        if (nproc.GT.1) return

        h = len_phys(1)/(snpart)
        ALLOCATE(p(ndim,npart))
        k = 0
        do i=1,snpart
            do j=1,snpart
                k = k + 1
                p(1,k) = h/2 + (i-1)*h
                p(2,k) = h/2 + (j-1)*h
            enddo
        enddo
        npart = k
        bcdef(1:2) = ppm_param_bcdef_periodic
        bcdef(3:4) = ppm_param_bcdef_symmetry
        bcdef(5:6) = ppm_param_bcdef_freespace
        !bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        CALL ppm_map_part_global(topoid,p,npart,info)
        CALL ppm_map_part_send(npart,newnpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart

        CALL ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)

        CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,1,gl,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)
        assert_equal(mpart-npart,104)
        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart,mpart)


    end test

    test symbc_and_map_load
        ! test ghost get map load and sym BC

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_util_dbg
          INTEGER                         :: npart = 2
          INTEGER                         :: newnpart
          INTEGER                         :: mpart
          REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
          REAL(MK),DIMENSION(:)  ,POINTER :: w => NULL()
          REAL(MK),DIMENSION(2)           :: check
          REAL(MK), PARAMETER             :: gl = 0.1_mk

        if (nproc.GT.1) return

        ALLOCATE(p(ndim,npart),w(npart))
        p(1,1) = 0.05_mk  ! left
        p(2,1) = 0.5_mk
        p(1,2) = 0.95_mk  ! right
        p(2,2) = 0.5_mk
        !p(1,3) = 0.5_mk   ! bottom
        !p(2,3) = 0.05_mk
        !p(1,4) = 0.5_mk   ! top
        !p(2,4) = 0.95_mk


        w(:) = 1.0_mk

        bcdef(1:2) = ppm_param_bcdef_symmetry
        bcdef(3:4) = ppm_param_bcdef_periodic

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        CALL ppm_map_part_global(topoid,p,npart,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,newnpart,info)
        CALL ppm_map_part_pop(w,npart,newnpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart
        CALL ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart)

        CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(w,npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,2,1,info,p,npart,mpart)
        CALL ppm_map_part_store(info)

        CALL ppm_map_part_ghost_put(topoid,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_ghost_pop(w,1,npart,mpart,info)

        p(1,1) = p(1,1)-0.01_mk
        p(2,1) = p(2,1)+0.1_mk
        p(1,2) = p(1,2)+0.01_mk
        p(2,2) = p(2,2)+0.1_mk

        CALL ppm_map_part_load(info)
        CALL ppm_map_part_push(p,ndim,npart,info,.TRUE.)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,3,1,info,p,npart,mpart)

        assert_equal(mpart-npart,2) ! check number of ghosts

        ! now go through all particles and try to find their ghosts

        check(1) = -0.04_mk  ! left
        check(2) =  0.6_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) =  1.04_mk  ! right
        check(2) =  0.6_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))


        DEALLOCATE(w)


    end test

    test symbc_and_map_load2
        ! another test ghost get map load and sym BC

        USE ppm_module_data
        USE ppm_module_topo_typedef
        USE ppm_module_mktopo
        USE ppm_module_topo_check
        USE ppm_module_util_dbg
          INTEGER                         :: npart = 2
          INTEGER                         :: newnpart
          INTEGER                         :: mpart
          REAL(MK),DIMENSION(:,:),POINTER :: p => NULL()
          REAL(MK),DIMENSION(:)  ,POINTER :: w => NULL()
          REAL(MK),DIMENSION(2)           :: check
          REAL(MK), PARAMETER             :: gl = 0.1_mk

        if (nproc.GT.1) return

        ALLOCATE(p(ndim,npart),w(npart))
        p(1,1) = 0.5_mk   ! bottom
        p(2,1) = 0.05_mk
        p(1,2) = 0.5_mk   ! top
        p(2,2) = 0.95_mk


        w(:) = 1.0_mk

        bcdef(1:2) = ppm_param_bcdef_freespace
        bcdef(3:4) = ppm_param_bcdef_symmetry

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        CALL ppm_map_part_global(topoid,p,npart,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,newnpart,info)
        CALL ppm_map_part_pop(w,npart,newnpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart
        CALL ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,1,1,info,p,npart)

        CALL ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(w,npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,2,1,info,p,npart,mpart)
        CALL ppm_map_part_store(info)

        CALL ppm_map_part_ghost_put(topoid,info)
        CALL ppm_map_part_push(w,npart,info)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_ghost_pop(w,1,npart,mpart,info)

        p(1,1) = p(1,1)+0.1_mk
        p(2,1) = p(2,1)+0.01_mk
        p(1,2) = p(1,2)-0.1_mk
        p(2,2) = p(2,2)+0.01_mk

        CALL ppm_map_part_load(info)
        CALL ppm_map_part_push(p,ndim,npart,info,.TRUE.)
        CALL ppm_map_part_send(npart,mpart,info)
        CALL ppm_map_part_pop(p,ndim,npart,mpart,info)

        CALL ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !CALL ppm_dbg_print(topoid,gl,3,1,info,p,npart,mpart)

        assert_equal(mpart-npart,2) ! check number of ghosts

        ! now go through all particles and try to find their ghosts

        check(1) =  0.6_mk  ! left
        check(2) = -0.06_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) =  0.4_mk  ! right
        check(2) =  1.04_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        DEALLOCATE(w)

    end test

end test_suite
