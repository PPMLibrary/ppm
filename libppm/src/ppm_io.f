      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_io
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs parallel I/O (read or write)
      !                 of scalars and 1d - 5d arrays on a distributed 
      !                 or centralized unit opened with ppm_io_open.
      !
      !  Input        : iUnit          (I) I/O unit (as returned by
      !                                    ppm_io_open).
      !                 actn           (I) I/O action. OPTIONAL. One of:
      !                                       ppm_param_io_read
      !                                       ppm_param_io_write
      !                                    If not specified, default is
      !                                    ppm_param_io_read
      !                 dist           (I) I/O distribution. OPTIONAL. 
      !                                    For read one of:
      !                                       ppm_param_io_same
      !                                       ppm_param_io_split
      !                                    If not specified, default is
      !                                    ppm_param_io_same
      !                                    For write one of:
      !                                       ppm_param_io_root
      !                                       ppm_param_io_concat
      !                                       ppm_param_io_sum
      !                                    If not specified, default is
      !                                    ppm_param_io_concat
      !                 iofmt          (C) Format string (for ASCII I/O).
      !                                    OPTIONAL. Default is (*)
      !                 ioprec         (I) Precision for binary I/O.
      !                                    One of:
      !                                       ppm_param_io_single
      !                                       ppm_param_io_double
      !                                    OPTIONAL. Default is
      !                                    leave-as-is.
      !
      !  Input/output : adata(...)     (O) data to be written/read.
      !                                    Overloaded types are: single,
      !                                    double, integer, single complex,
      !                                    double complex, logical,
      !                                    character. Can be either scalar,
      !                                    1d, 2d, 3d, 4d or 5d. Character
      !                                    strings can only be scalar (i.e.
      !                                    no array of strings).
      !
      !  Output       : stat           (I) return code. 0 on success.
      !                                    OPTIONAL
      !
      !  Remarks      : Many arguments are OPTIONAL. It is adivised
      !                 to use explicit argument naming when using this
      !                 routine. Example:
      !                 CALL ppm_io(20,xp,ACTN=ppm_param_io_read,STAT=info)
      !
      !                 Beware: RESHAPE operates on the stack, so on machines
      !                 with limited stack size this limits the max possible
      !                 data size. Therefore we use DO loops (not slower) to
      !                 remain on the heap. On vector machines we still use
      !                 RESHAPE because it is faster.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/28 07:17:13  ivos
      !  Added explicit DO loops as an alternative to RESHAPE (see Remarks).
      !
      !  Revision 1.12  2006/09/04 18:34:49  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.10  2004/11/11 15:24:31  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.9  2004/10/01 16:09:00  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.8  2004/07/26 07:46:36  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.7  2004/07/16 14:46:25  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.6  2004/05/27 14:49:58  ivos
      !  Bugfix: info was not caught after CALL to ppm_doio.
      !
      !  Revision 1.5  2004/05/27 10:41:56  ivos
      !  Complete rewrite: now converts everything to a 1d array and the whole
      !  MPI and I/O stuff is outsourced to ppm_doio.
      !
      !  Revision 1.4  2004/05/18 09:05:22  ivos
      !  Serial version did not compile. fixed.
      !
      !  Revision 1.3  2004/05/14 14:50:03  ivos
      !  bugfix: changed WRITEs of CHARACTERs to use TRIM in order to avoid
      !  chunk being written after the end of the string. Also fixed wrong
      !  I/O format strings for list directed I/O.
      !
      !  Revision 1.2  2004/05/13 11:39:09  ivos
      !  Added proper error handling for the parallel cases as well as proper
      !  handling of character strings using MPI.
      !
      !  Revision 1.1  2004/05/06 10:44:59  ivos
      !  Preliminary check-in. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __DIM == 0
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_0ds(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_0dd(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_0dsc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_0ddc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_0di(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_0dl(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __CHARACTER
      SUBROUTINE ppm_io_0dc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#endif

#elif __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_1ds(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_1dd(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_1dsc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_1ddc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_1di(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_1dl(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#endif

#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_2ds(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_2dd(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_2dsc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_2ddc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_2di(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_2dl(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#endif

#elif __DIM == 3
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_3ds(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_3dd(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_3dsc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_3ddc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_3di(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_3dl(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#endif

#elif __DIM == 4
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_4ds(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_4dd(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_4dsc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_4ddc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_4di(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_4dl(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#endif

#elif __DIM == 5
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_5ds(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_5dd(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_5dsc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_5ddc(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_5di(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_5dl(iUnit,adata,actn,dist,iofmt,ioprec,stat)
#endif
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_doio
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                              , INTENT(IN   ) :: iUnit
#if   __DIM == 0
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                             , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)                          , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER                              , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL                              , INTENT(INOUT) :: adata
#elif __KIND == __CHARACTER
      CHARACTER(LEN=*)                     , INTENT(INOUT) :: adata
#endif

#elif __DIM == 1
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)        , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(INOUT) :: adata
#endif

#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:,:)        , INTENT(INOUT) :: adata
#endif

#elif __DIM == 3
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:,:,:)        , INTENT(INOUT) :: adata
#endif

#elif __DIM == 4
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:,:,:,:)        , INTENT(INOUT) :: adata
#endif

#elif __DIM == 5
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:,:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:,:,:,:)        , INTENT(INOUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:,:,:,:,:)        , INTENT(INOUT) :: adata
#endif
#endif
      INTEGER         , OPTIONAL, INTENT(IN   ) :: actn,dist,ioprec
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: iofmt
      INTEGER         , OPTIONAL, INTENT(  OUT) :: stat
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0
      INTEGER                          :: info,iactn,idist,iprec
      INTEGER                          :: i,j,k,l,m
      CHARACTER(LEN=9)                 :: subname
      CHARACTER(LEN=ppm_char)          :: ifmt
      INTEGER, DIMENSION(5)            :: lda
      INTEGER, DIMENSION(1)            :: ldc
      INTEGER                          :: iopt,ndata,ishift,idata
      LOGICAL                          :: lchar
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:), POINTER :: abuf
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:), POINTER :: abuf
#elif __KIND == __INTEGER | __KIND == __CHARACTER
      INTEGER    , DIMENSION(:), POINTER :: abuf
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:), POINTER :: abuf
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      WRITE(subname,'(A,I1.1,A)') 'ppm_io_',__DIM,'d'
      CALL substart(subname,t0,info)
      lchar = .FALSE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check and decode arguments
      !-------------------------------------------------------------------------
      IF (iUnit .LE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,subname,    &
     &        'Unit number needs to be > 0',__LINE__,info)
          GOTO 9999
      ENDIF

      iactn = ppm_param_io_read
      IF (PRESENT(actn)) THEN
          IF     ((actn .NE. ppm_param_io_read) .AND.    &
     &            (actn .NE. ppm_param_io_write)) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,subname,    &
     &           'Invalid I/O action specified.',__LINE__,info)
             GOTO 9999
          ENDIF
          iactn = actn
      ENDIF

      IF (iactn .EQ. ppm_param_io_read) THEN
          idist = ppm_param_io_same
          IF (PRESENT(dist)) THEN
              IF     ((dist .NE. ppm_param_io_same) .AND.    &
     &                (dist .NE. ppm_param_io_split)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,subname,    &
     &               'Invalid I/O distribution specified.',__LINE__,info)
                 GOTO 9999
              ENDIF
              idist = dist
          ENDIF
      ELSE
          idist = ppm_param_io_concat
          IF (PRESENT(dist)) THEN
              IF     ((dist .NE. ppm_param_io_root) .AND.    &
     &                (dist .NE. ppm_param_io_concat) .AND.  &
     &                (dist .NE. ppm_param_io_sum)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,subname,    &
     &               'Invalid I/O distribution specified.',__LINE__,info)
                 GOTO 9999
              ENDIF
              idist = dist
          ENDIF
      ENDIF

#if __KIND == __CHARACTER
      IF (idist .EQ. ppm_param_io_split) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,subname,    &
     &     'Split distribution not possible for characters.',__LINE__,info)
          GOTO 9999 
      ENDIF
      IF (idist .EQ. ppm_param_io_sum) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,subname,    &
     &     'Sum reduction not possible for characters.',__LINE__,info)
          GOTO 9999 
      ENDIF
#endif

      ifmt = ''
      IF (PRESENT(iofmt)) THEN
          i = LEN_TRIM(iofmt)
          ifmt(1:i) = iofmt(1:i)
      ENDIF

      iprec = ppm_param_undefined
      IF (PRESENT(ioprec)) THEN
          IF     ((ioprec .NE. ppm_param_io_single) .AND.    &
     &            (ioprec .NE. ppm_param_io_double)) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,subname,    &
     &           'Invalid I/O precision specified.',__LINE__,info)
             GOTO 9999
          ENDIF
          iprec = ioprec
      ENDIF

      !-------------------------------------------------------------------------
      !  Point to proper work storage
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION 
      abuf => abuf_s
#elif __KIND == __DOUBLE_PRECISION
      abuf => abuf_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX 
      abuf => abuf_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      abuf => abuf_dc
#elif __KIND == __INTEGER | __KIND == __CHARACTER
      abuf => abuf_i
#elif __KIND == __LOGICAL
      abuf => abuf_l
#endif

#if __KIND == __CHARACTER
      !-------------------------------------------------------------------------
      !  Convert character string to integer vector
      !-------------------------------------------------------------------------
      ndata = LEN(adata)

      iopt = ppm_param_alloc_fit
      ldc(1) = ndata
      CALL ppm_alloc(abuf,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,subname,'1d buffer ABUF',__LINE__,info)
          GOTO 9999
      ENDIF

      IF (iactn .EQ. ppm_param_io_write) THEN
          ishift = IACHAR('a')
          ! convert to integers relative to the lower-case a (in order to
          ! be independent of ASCII table shifts on the different
          ! processors)
          DO i=1,ndata
              abuf(i) = IACHAR(adata(i:i))-ishift
          ENDDO
      ENDIF

      lchar = .TRUE.

#elif __DIM == 0
      !-------------------------------------------------------------------------
      !  Convert to 1d array
      !-------------------------------------------------------------------------
      ndata = 1

      iopt = ppm_param_alloc_fit
      ldc(1) = ndata
      CALL ppm_alloc(abuf,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,subname,'1d buffer ABUF',__LINE__,info)
          GOTO 9999
      ENDIF

      IF (iactn .EQ. ppm_param_io_write) THEN
          abuf(1) = adata
      ENDIF

#elif __DIM > 1
      !-------------------------------------------------------------------------
      !  Convert to 1d array
      !-------------------------------------------------------------------------
      ndata = 1
      DO i=1,__DIM
          lda(i) = SIZE(adata,i)
          ndata = ndata*lda(i)
      ENDDO

      iopt = ppm_param_alloc_fit
      ldc(1) = ndata
      CALL ppm_alloc(abuf,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,subname,'1d buffer ABUF',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (iactn .EQ. ppm_param_io_write) THEN
#ifdef __VECTOR
          !---------------------------------------------------------------------
          !  On vector machines RESHAPE is the fastest
          !---------------------------------------------------------------------
          abuf = RESHAPE(adata,ldc)
#else

          !---------------------------------------------------------------------
          !  On scalar machines use DO loops since RESHAPE operates on the 
          !  stack and we easily run into the limit. DO loops avoid this.
          !---------------------------------------------------------------------
#if   __DIM == 2
          idata = 0
          DO i=1,lda(2)
              DO j=1,lda(1)
                  idata = idata + 1
                  abuf(idata) = adata(j,i)
              ENDDO
          ENDDO
#elif __DIM == 3
          idata = 0
          DO i=1,lda(3)
              DO j=1,lda(2)
                  DO k=1,lda(1)
                      idata = idata + 1
                      abuf(idata) = adata(k,j,i)
                  ENDDO
              ENDDO
          ENDDO
#elif __DIM == 4
          idata = 0
          DO i=1,lda(4)
              DO j=1,lda(3)
                  DO k=1,lda(2)
                      DO l=1,lda(1)
                          idata = idata + 1
                          abuf(idata) = adata(l,k,j,i)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
#elif __DIM == 5
          idata = 0
          DO i=1,lda(5)
              DO j=1,lda(4)
                  DO k=1,lda(3)
                      DO l=1,lda(2)
                          DO m=1,lda(1)
                              idata = idata + 1
                              abuf(idata) = adata(m,l,k,j,i)
                          ENDDO
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
#endif
#endif
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Do the I/O
      !-------------------------------------------------------------------------
#if __DIM == 1
      CALL ppm_doio(iUnit,adata,iactn,idist,ifmt,iprec,lchar,__DIM,info)
#else
      CALL ppm_doio(iUnit,abuf,iactn,idist,ifmt,iprec,lchar,__DIM,info)
#endif
      IF (info .NE. 0) GOTO 9999

#if __KIND == __CHARACTER
      !-------------------------------------------------------------------------
      !  Convert to character string
      !-------------------------------------------------------------------------
      IF (iactn .EQ. ppm_param_io_read) THEN
          ishift = IACHAR('a') 
          DO i=1,ndata
              adata(i:i) = CHAR(abuf(i)+ishift)
          ENDDO
      ENDIF

      iopt = ppm_param_dealloc
      CALL ppm_alloc(abuf,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,subname,'1d buffer ABUF',__LINE__,info)
      ENDIF

#elif __DIM == 0
      !-------------------------------------------------------------------------
      !  Convert to scalar
      !-------------------------------------------------------------------------
      IF (iactn .EQ. ppm_param_io_read) THEN
          adata = abuf(1)
      ENDIF

      iopt = ppm_param_dealloc
      CALL ppm_alloc(abuf,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,subname,'1d buffer ABUF',__LINE__,info)
      ENDIF

#elif __DIM > 1
      !-------------------------------------------------------------------------
      !  Convert back to array
      !-------------------------------------------------------------------------
      IF (iactn .EQ. ppm_param_io_read) THEN
#ifdef __VECTOR
          !---------------------------------------------------------------------
          !  On vector machines RESHAPE is the fastest
          !---------------------------------------------------------------------
          adata = RESHAPE(abuf,lda(1:__DIM))
#else

          !---------------------------------------------------------------------
          !  On scalar machines use DO loops since RESHAPE operates on the 
          !  stack and we easily run into the limit. DO loops avoid this.
          !---------------------------------------------------------------------
#if   __DIM == 2
          idata = 0
          DO i=1,lda(2)
              DO j=1,lda(1)
                  idata = idata + 1
                  adata(j,i) = abuf(idata)
              ENDDO
          ENDDO
#elif __DIM == 3
          idata = 0
          DO i=1,lda(3)
              DO j=1,lda(2)
                  DO k=1,lda(1)
                      idata = idata + 1
                      adata(k,j,i) = abuf(idata)
                  ENDDO
              ENDDO
          ENDDO
#elif __DIM == 4
          idata = 0
          DO i=1,lda(4)
              DO j=1,lda(3)
                  DO k=1,lda(2)
                      DO l=1,lda(1)
                          idata = idata + 1
                          adata(l,k,j,i) = abuf(idata)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
#elif __DIM == 5
          idata = 0
          DO i=1,lda(5)
              DO j=1,lda(4)
                  DO k=1,lda(3)
                      DO l=1,lda(2)
                          DO m=1,lda(1)
                              idata = idata + 1
                              adata(m,l,k,j,i) = abuf(idata)
                          ENDDO
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
#endif
#endif
      ENDIF

      iopt = ppm_param_dealloc
      CALL ppm_alloc(abuf,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,subname,'1d buffer ABUF',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      IF (PRESENT(stat)) stat = info
      CALL substop(subname,t0,info)
      RETURN
#if   __DIM == 0
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_0ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_0dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_0dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_0ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_0di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_0dl
#elif __KIND == __CHARACTER
      END SUBROUTINE ppm_io_0dc
#endif

#elif __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_1ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_1ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_1di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_1dl
#endif

#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_2dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_2ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_2di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_2dl
#endif

#elif __DIM == 3
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_3dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_3dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_3ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_3di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_3dl
#endif

#elif __DIM == 4
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_4ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_4dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_4dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_4ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_4di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_4dl
#endif

#elif __DIM == 5
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_5ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_5dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_5dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_5ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_5di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_5dl
#endif
#endif
