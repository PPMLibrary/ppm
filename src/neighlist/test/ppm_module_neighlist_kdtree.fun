test_suite ppm_module_neighlist_kdtree

#ifdef __MPI
  INCLUDE "mpif.h"
#endif

  INTEGER,  PARAMETER                 :: debug = 0
  INTEGER,  PARAMETER                 :: mk = kind(1.0d0) !kind(1.0e0)
  REAL(MK), PARAMETER                 :: skin = 0._mk
  INTEGER,  PARAMETER                 :: ndim=3
  INTEGER                             :: decomp,assig,tolexp
  REAL(MK)                            :: tol,min_rcp,max_rcp
  INTEGER                             :: info,comm,rank,nproc
  INTEGER                             :: topoid
  INTEGER                             :: np = 100000
  INTEGER                             :: npart
  INTEGER                             :: mp
  INTEGER                             :: newnp
  REAL(MK), DIMENSION(:,:), POINTER   :: xp => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: rcp => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: min_phys => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: max_phys => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: h => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: len_phys => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: ghostlayer => NULL()
  INTEGER,  DIMENSION(:  ), POINTER   :: ghostsize => NULL()
  INTEGER                             :: i,ip,ineigh
  INTEGER,  DIMENSION(6)              :: bcdef
  REAL(MK), DIMENSION(:  ), POINTER   :: cost => NULL()
  INTEGER,  DIMENSION(:  ), POINTER   :: nm => NULL()
  INTEGER                             :: seedsize
  INTEGER,  DIMENSION(:), ALLOCATABLE :: seed
  REAL(MK), DIMENSION(:), ALLOCATABLE :: randnb
  INTEGER                             :: isymm = 0
  LOGICAL                             :: lsymm = .false.,ok
  REAL(MK)                            :: t0,t1,t2,t3
  REAL(MK)                            :: eps

  init
    USE ppm_module_data
    USE ppm_module_topo_typedef
    USE ppm_module_init

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
    &       ghostsize(ndim),ghostlayer(2*ndim),&
    &       nm(ndim),h(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_mk
    max_phys(1:ndim) = 1.0_mk
    len_phys(1:ndim) = max_phys-min_phys
    ghostsize(1:ndim) = 2
    bcdef(1:6) = ppm_param_bcdef_periodic

    eps = EPSILON(1.0_MK)
    tolexp = INT(LOG10(EPSILON(1.0_MK)))

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

    DEALLOCATE(min_phys,max_phys,len_phys,ghostsize,nm,h)
  end finalize


  setup
    CALL RANDOM_SEED(SIZE=seedsize)
    ALLOCATE(seed(seedsize))
    ALLOCATE(randnb((1+ndim)*np),STAT=info)
    DO i=1,seedsize
       seed(i)=10+i*i*(rank+1)
    ENDDO
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(randnb)

    ALLOCATE(xp(ndim,np),rcp(np),STAT=info)
  end setup

  teardown
    DEALLOCATE(xp,rcp,STAT=info)
    DEALLOCATE(seed,randnb)
  end teardown

  test kdtree_test
    USE ppm_module_kdtree
    USE ppm_module_inl_k_vlist
    USE ppm_module_mktopo
    USE ppm_module_particles_typedef
    USE ppm_module_util_time

    REAL(MK), DIMENSION(ndim) :: coord

    TYPE(kdtree2_d), POINTER :: tree
    TYPE(kdtree2_result_d), DIMENSION(:), ALLOCATABLE, TARGET :: results,resultsb
    TYPE(ppm_t_particles_d), POINTER :: Part
    CLASS(ppm_t_neighlist_d_), POINTER :: Nlist

    REAL(MK) :: cutoff,maxdeviation

    INTEGER :: knn

    start_subroutine("kdtree_test")

    !----------------
    ! make topology
    !----------------
    decomp =  ppm_param_decomp_cuboid
    assig  = ppm_param_assign_internal

    topoid = 0
    cutoff = 0.15_mk
    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
    Assert_Equal(info,0)

    ALLOCATE(Part,STAT=info)
    Assert_Equal(info,0)

    CALL Part%initialize(np,info,topoid=topoid, &
    &    distrib=ppm_param_part_init_random,name="Part")
    Assert_Equal(info,0)

    Part%ghostlayer=.06_MK

    CALL Part%map(info,global=.TRUE.,topoid=topoid)
    Assert_Equal(info,0)

    CALL Part%map_ghosts(info)
    Assert_Equal(info,0)

    ALLOCATE(tree,STAT=info)
    Assert_Equal(info,0)

    CALL Part%get_xp(xp,info,with_ghosts=.TRUE.)
    Assert_Equal(info,0)

    CALL ppm_util_time(t0)
    CALL ppm_util_time(t0)
    ! this is how you create a tree.
    CALL tree%create(xp,info,dim=ndim,sort=.FALSE.,rearrange=.FALSE.)
    Assert_Equal(info,0)
    CALL ppm_util_time(t1)
    stdout('REAL(Part%Mpart,MK)/(t1-t0)'," points per second built for non-rearranged tree.")
    CALL tree%destroy(info)
    Assert_Equal(info,0)

    CALL ppm_util_time(t0)
    ! this is how you create a tree.
    CALL tree%create(xp,info,dim=ndim,sort=.FALSE.,rearrange=.TRUE.)
    Assert_Equal(info,0)
    CALL ppm_util_time(t1)
    stdout('REAL(Part%Mpart,MK)/(t1-t0)'," points per second built for rearranged tree.")
    CALL tree%destroy(info)
    Assert_Equal(info,0)

    CALL ppm_util_time(t0)
    ! this is how you create a tree.
    CALL tree%create(xp,info,dim=ndim,sort=.TRUE.,rearrange=.TRUE.)
    Assert_Equal(info,0)
    CALL ppm_util_time(t1)
    stdout('REAL(Part%Mpart,MK)/(t1-t0)'," points per second built for sorted rearranged tree.")
    CALL tree%destroy(info)
    Assert_Equal(info,0)



    CALL ppm_util_time(t0)
    CALL Part%create_neighlist(Part,info,skin=0.0_MK)
    CALL ppm_util_time(t1)
    CALL Part%comp_neighlist(info,incl_ghosts=.FALSE.)
    CALL ppm_util_time(t2)
    stdout('t1-t0'," seconds to create neighlist.")
    stdout('t2-t1'," seconds to compute Verlet list.")
    Assert_Equal(info,0)

    Assert_true(Part%has_neighlist())

    Nlist => Part%get_neighlist()

    knn=50
    ALLOCATE(results(knn+1),resultsb(knn+1),STAT=info)
    Assert_Equal(info,0)

    CALL ppm_util_time(t0)
    CALL tree%create(xp,info,dim=ndim,sort=.TRUE.,rearrange=.TRUE.)
    CALL ppm_util_time(t1)
    DO ip=1,Part%Npart
       ! If the tree is not sorted you need to remove the self
       ! particle from list of neighbors
       CALL kdtree2_n_nearest(tree,xp(1:ndim,ip), &
       &    knn+1,results,info)
       Assert_Equal(info,0)
       Nlist%vlist(1:knn,ip)=results(2:knn+1)%idx
    ENDDO
    Nlist%nvlist=knn
    CALL ppm_util_time(t2)
    stdout('t1-t0'," seconds to create the tree.")
    stdout('t2-t1'," seconds to create Verlet list using the tree.")
    Assert_Equal(info,0)

    !compare the results of kdtree nearest search with brute search
    DO ip=1,10  !Part%Npart
       CALL kdtree2_n_nearest(tree,xp(1:ndim,ip), &
       &    knn+1,results,info)
       Assert_Equal(info,0)

       CALL kdtree2_n_nearest_brute_force(tree,xp(1:ndim,ip), &
       &    knn+1,resultsb,info)
       Assert_Equal(info,0)

       maxdeviation=MAXVAL(ABS(results(1:knn+1)%dis-resultsb(1:knn+1)%dis))
       IF (ANY(results(1:knn+1)%idx.NE.resultsb(1:knn+1)%idx).OR. &
          (maxdeviation.GT.1.0E-8_MK)) THEN
          stdout("MISMATCH! @ particle ",ip)
       ENDIF
    ENDDO

    Nlist => NULL()
    CALL Part%set_xp(xp,info)
    Assert_Equal(info,0)

    DEALLOCATE(results,STAT=info)
    Assert_Equal(info,0)

    CALL tree%destroy(info)
    Assert_Equal(info,0)

    CALL Part%destroy(info)
    Assert_Equal(info,0)

    end_subroutine()
  end test

end test_suite