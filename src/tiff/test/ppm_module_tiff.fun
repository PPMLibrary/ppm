test_suite ppm_module_tiff

#ifdef __MPI
  INCLUDE "mpif.h"
#endif

  INTEGER,  PARAMETER                 :: debug = 0
  INTEGER,  PARAMETER                 :: MK = kind(1.0D0)
  INTEGER,  PARAMETER                 :: ndim=2
  INTEGER,  DIMENSION(4)              :: bcdef
  INTEGER                             :: decomp,assig,tolexp
  INTEGER                             :: info,comm,rank,nproc
  INTEGER                             :: topoid
  INTEGER, DIMENSION(3)               :: Ngrid
  REAL(MK), DIMENSION(:  ), POINTER   :: min_phys => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: max_phys => NULL()
  REAL(MK), DIMENSION(:  ), POINTER   :: cost => NULL()
  REAL(MK), DIMENSION(2)              :: Offset
  REAL(MK), DIMENSION(:,:), POINTER   :: wp => NULL()
  REAL(MK), DIMENSION(:),   POINTER   :: wpt => NULL()
  INTEGER,  DIMENSION(:  ), POINTER   :: nm => NULL()
  INTEGER,  DIMENSION(:  ), POINTER   :: ghostsize => NULL()
  INTEGER                             :: i,j,istrip

  init
    USE ppm_module_data
    USE ppm_module_init

#ifdef __MPI
    comm = MPI_COMM_WORLD
    CALL MPI_Comm_rank(comm,rank,info)
    CALL MPI_Comm_size(comm,nproc,info)
#else
    rank = 0
    nproc = 1
#endif

    bcdef = ppm_param_bcdef_freespace
    tolexp = INT(LOG10(EPSILON(1._mk)))+10

    CALL ppm_init(ndim,MK,tolexp,comm,debug,info,99)

    ALLOCATE(min_phys(ndim),max_phys(ndim),ghostsize(ndim),STAT=info)

    min_phys(1:ndim) = 0.0_MK
    ghostsize(1:ndim)= 0
    Offset=0.0_MK
  end init

  finalize
    USE ppm_module_finalize
    CALL ppm_finalize(info)
    DEALLOCATE(min_phys,max_phys,ghostsize)
  end finalize

  setup
  end setup

  teardown
  end teardown

  test tiff
    ! test tiff module

    USE ppm_module_data
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_write
    USE ppm_module_topo_typedef
    USE ppm_module_mktopo
    USE ppm_module_tiff
    USE ppm_module_interfaces
    USE ppm_module_mesh_typedef
    USE ppm_module_field_typedef

    CLASS(ppm_t_equi_mesh_), POINTER :: mesh
    CLASS(ppm_t_field_),     POINTER :: image
    TYPE(ppm_t_topo),        POINTER :: topo
    CLASS(ppm_t_subpatch_),  POINTER :: sbpitr

#ifdef __MPI3
    INTEGER :: request
#endif

    CHARACTER(LEN=ppm_char) :: inputimage
    CHARACTER(LEN=ppm_char) :: outputimage

    start_subroutine("tiff")

#ifdef __TIFF
    WRITE(inputimage,'(A,A1)') "src/tiff/test/icecream.tif",CHAR(0)

    !Reading TIFF header file to know the image size
    !This part will broadcast the header information to all processors
    info=read_tiff_info(inputimage,Ngrid)
    Assert_Equal(info,0)

    IF (rank.EQ.0) THEN
       stdout("width=",'Ngrid(1)',"length=",'Ngrid(2)')
       stdout("bitsPerSample=",bitsPerSample,"samplesPerPixel=",samplesPerPixel)
    ENDIF

    !maximum physical size
    max_phys(1:ndim) = REAL(Ngrid(1:ndim)-1,MK)+EPSILON(1.0_MK)

    !----------------
    ! make topology
    !----------------
    assig  = ppm_param_assign_internal
    decomp = ppm_param_decomp_xpencil

    topoid = 0

    CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,&
    &    bcdef,REAL(MAXVAL(ghostsize),MK),cost,info)
    Assert_Equal(info,0)

    !----------------
    ! make mesh to read image file
    !----------------
    ALLOCATE(ppm_t_equi_mesh::mesh,STAT=info)
    Assert_Equal(info,0)

    CALL mesh%create(topoid,Offset,info, &
    &    Nm=Ngrid,ghostsize=ghostsize)
    Assert_Equal(info,0)

    CALL mesh%def_uniform(info)
    Assert_Equal(info,0)

    ALLOCATE(ppm_t_field::image,STAT=info)
    Assert_Equal(info,0)

    CALL image%create(1,info,name="image_pixel_intensity")
    Assert_Equal(info,0)

    CALL image%discretize_on(mesh,info)
    Assert_Equal(info,0)

    info=open_tiff(inputimage)
    Assert_Equal(info,0)

    sbpitr => mesh%subpatch%begin()
    DO WHILE (ASSOCIATED(sbpitr))
       nm => sbpitr%nnodes

       CALL sbpitr%get_field(image,wp,info)
       Assert_Equal(info,0)

       j=1
       DO istrip=sbpitr%istart(ndim),sbpitr%iend(ndim)
          !reading part of the tiff file associated with this processor line by line
          info=read_tiff_scanline(wp(1:nm(1),j),istrip-1,bitsPerSample,0,nm(1))
          Assert_Equal(info,0)

          j=j+1
       ENDDO

       sbpitr => mesh%subpatch%next()
    ENDDO

    info=close_tiff()
    Assert_Equal(info,0)

    sbpitr => mesh%subpatch%begin()
    DO WHILE (ASSOCIATED(sbpitr))
       nm => sbpitr%nnodes

       CALL sbpitr%get_field(image,wp,info)
       Assert_Equal(info,0)

       WRITE(outputimage,'(A,I0,A1,I0,A4,A1)') "src/tiff/test/IO_",sbpitr%isub,"_",rank,".tif",CHAR(0)

       OPEN(UNIT=10000,FILE=outputimage,IOSTAT=info)
       IF (info.EQ.0) CLOSE(10000,STATUS='DELETE')

       info=open_write_tiff(outputimage)
       Assert_Equal(info,0)

       info=write_tiff_header(bitsPerSampleW,1,nm(1),nm(2))
       Assert_Equal(info,0)

       !Write the output to an image line by line
       !DO j=1,nm(ndim)
       !   info=write_tiff_scanline(wp(1:nm(1),j),j-1,bitsPerSampleW,0,nm(1))
       !   Assert_Equal(info,0)
       !ENDDO

       !Write the strip of data to output
       info=write_tiff_strip(wp(1:nm(1),1:nm(2)),bitsPerSampleW,nm(1),nm(2),0,1)
       Assert_Equal(info,0)

       OPEN(UNIT=10000,FILE=outputimage,IOSTAT=info)
       IF (info.EQ.0) CLOSE(10000,STATUS='DELETE')

       info=close_tiff()
       Assert_Equal(info,0)

       sbpitr => mesh%subpatch%next()
    ENDDO

    CALL mesh%destroy(info)
    Assert_Equal(info,0)

    dealloc_pointer("mesh")

#else
    Assert_Equal(info,0)
#endif

    end_subroutine()
  end test

end test_suite
