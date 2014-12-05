test_suite ppm_module_vbp_particles

  USE ppm_module_particles_typedef
  USE ppm_module_vbp_typedef
  USE ppm_module_topo_typedef
  USE ppm_module_field_typedef
  USE ppm_module_operator_typedef
  USE ppm_module_interfaces
  USE ppm_module_data

#ifdef __MPI
  INCLUDE "mpif.h"
#endif

  INTEGER, PARAMETER              :: debug = 0
  INTEGER, PARAMETER              :: MK = ppm_kind_double
  REAL(MK),PARAMETER              :: tol=EPSILON(1._MK)*100
  REAL(MK),PARAMETER              :: pi = ACOS(-1._MK)
  REAL(MK),PARAMETER              :: skin = 0._MK
  INTEGER,PARAMETER               :: ndim=2
  INTEGER                         :: decomp,assig,tolexp
  INTEGER                         :: info,comm,rank,nproc,topoid
  INTEGER                         :: np_global
  REAL(MK), DIMENSION(:),   POINTER :: cutoff
  REAL(MK), DIMENSION(:,:), POINTER :: xp => NULL()
  REAL(MK), DIMENSION(:  ), POINTER :: min_phys, max_phys
  REAL(MK), DIMENSION(:  ), POINTER :: len_phys
  INTEGER                         :: i,j,k,ip,wp_id,iopt
  INTEGER                         :: nstep
  INTEGER, DIMENSION(3)           :: ldu
  INTEGER, DIMENSION(2*ndim)      :: bcdef
  REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
  INTEGER                         :: isymm = 0
  REAL(MK)                        :: t0,t1,t2,t3
  INTEGER                         :: seedsize
  INTEGER,  DIMENSION(:), ALLOCATABLE :: seed
  INTEGER, DIMENSION(:),   POINTER :: nvlist => NULL()
  INTEGER, DIMENSION(:,:), POINTER :: vlist  => NULL()
  REAL(MK)                        :: err

  init
    USE ppm_module_init
    USE ppm_module_mktopo

    ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    max_phys(1:ndim) = (/1.01_MK,0.65_MK/)
    len_phys(1:ndim) = max_phys-min_phys
    bcdef = ppm_param_bcdef_freespace

#ifdef __MPI
    CALL MPI_Comm_dup(MPI_COMM_WORLD,comm,info)
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    rank = 0
    nproc = 1
#endif
    tolexp = INT(LOG10(EPSILON(1._MK)))+10
    CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

    !----------------
    ! make topology
    !----------------
    decomp=ppm_param_decomp_cuboid
    assig =ppm_param_assign_internal

    topoid = 0
    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,0.1_MK,cost,info)
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

  test neighlists
    USE ppm_module_util_time
    USE ppm_module_io_vtk

    TYPE(ppm_t_vbp_d), TARGET :: Part1

    CLASS(ppm_t_neighlist_d_), POINTER :: Nlist => NULL()
    CLASS(ppm_t_discr_data),   POINTER :: prop => NULL()

    REAL(MK) :: t0,t1,t2


    start_subroutine("neighlists")

    !--------------------------
    ! Create particle set and read the x&y coordinates and cutoffs
    ! from the input file
    !--------------------------

    IF (rank.EQ.0) THEN
       np_global = 61649

       CALL Part1%create(np_global,info,name="Part1")
       Assert_Equal(info,0)

       iopt=ppm_param_alloc_fit
       ldu(1)=np_global
       CALL ppm_alloc(cutoff,ldu,iopt,info)
       Assert_Equal(info,0)

       OPEN(UNIT=100,FILE="src/vbp/test/x.txt",STATUS='OLD')
       OPEN(UNIT=200,FILE="src/vbp/test/y.txt",STATUS='OLD')
       OPEN(UNIT=300,FILE="src/vbp/test/r.txt",STATUS='OLD')

       CALL Part1%get_xp(xp,info)
       Assert_Equal(info,0)

       DO i=1,np_global
          READ(UNIT=100,FMT=*) xp(1,i)
          READ(UNIT=200,FMT=*) xp(2,i)
          READ(UNIT=300,FMT=*) cutoff(i)
       ENDDO

       CALL Part1%set_xp(xp,info)
       Assert_Equal(info,0)

       CLOSE(UNIT=100)
       CLOSE(UNIT=200)
       CLOSE(UNIT=300)
    ELSE
       np_global = 0

       CALL Part1%create(np_global,info,name="Part1")
       Assert_Equal(info,0)

       iopt=ppm_param_alloc_fit
       ldu(1)=np_global
       CALL ppm_alloc(cutoff,ldu,iopt,info)
       Assert_Equal(info,0)
    ENDIF

    CALL Part1%create_neighlist(Part1,info,skin=0.0_MK)
    Assert_Equal(info,0)

    CALL Part1%set_varying_cutoff(cutoff,info)
    Assert_Equal(info,0)

    iopt=ppm_param_dealloc
    CALL ppm_alloc(cutoff,ldu,iopt,info)
    Assert_Equal(info,0)

    !The hack to make it work
    !TODO
    !find a better way to avoid this hack
    Part1%active_topoid=topoid

    CALL Part1%apply_bc(info)
    Assert_Equal(info,0)

    global_mapping(Part1, topoid)

    CALL Part1%get(Part1%rcp,cutoff,info)
    Assert_Equal(info,0)

    Part1%ghostlayer = MAXVAL(cutoff(1:Part1%Npart))

    CALL Part1%set(Part1%rcp,cutoff,info)
    Assert_Equal(info,0)

    CALL Part1%map_ghosts(info)
    Assert_Equal(info,0)

    CALL ppm_util_time(t1)
    CALL Part1%comp_neighlist(info,incl_ghosts=.FALSE.)
    CALL ppm_util_time(t2)
    Assert_Equal(info,0)

    stdout("It took: ",'t2-t1'," secs")

    ! print particles to a VTK file
    !CALL ppm_vtk_particles("output",Part1,info)
    !Assert_Equal(info,0)

    CALL Part1%destroy(info)
    Assert_Equal(info,0)

    end_subroutine()
  end test

end test_suite
