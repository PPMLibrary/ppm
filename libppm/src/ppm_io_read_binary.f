      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_io_read_binary
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine provides the atomic reading
      !                 functionality for binary 1d arrays.
      !                 It is never used by the user program directly.
      !
      !  Input        : iUnit          (I) Fortran I/O unit to be used
      !                 nrec           (I) length of record to be read
      !                 ioprec         (I) Precision for I/O.
      !                                    One of:
      !                                       ppm_param_io_single
      !                                       ppm_param_io_double
      !                                    Integer and Logical cases are
      !                                    not influenced by this. No
      !                                    conversions are done if 
      !                                    any value other than the above
      !                                    is given.
      !                 lchar          (L) .TRUE. if a CHARACTER string is
      !                                    to be read and has to be
      !                                    converted to an INTEGER array
      !                                    before returning it.
      !
      !  Input/output : 
      !
      !  Output       : adata          (O) Data vector to be read. 1d
      !                                    array of either single, double,
      !                                    single complex, double complex,
      !                                    integer or logical elements.
      !                                    Needs to be allocated to proper
      !                                    size before calling this
      !                                    routine. Its size determines the
      !                                    number of elements read from the
      !                                    file.
      !                 info           (I) return code. 0 on success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_read_binary.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2006/02/03 09:35:19  ivos
      !  Fixed bug 00004: cbuf was only allocated in the __INTEGER case, but
      !  used in an else-branch. Added separate ifdefs for the __INTEGER
      !  case to make sure cbuf is ONLY used there.
      !
      !  Revision 1.6  2005/01/05 16:25:46  ivos
      !  bugfix: precision buffers were never deallocated.
      !
      !  Revision 1.5  2004/10/01 16:09:02  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/09/07 14:29:25  ivos
      !  fix: field length descriptors inserted in debug output in order to
      !  avoid compiler warnings on XLF.
      !
      !  Revision 1.3  2004/07/26 07:45:27  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.2  2004/05/27 14:52:10  ivos
      !  Bugfix: fixed misplaced GOTO
      !
      !  Revision 1.1  2004/05/27 10:43:42  ivos
      !  Initial implementation. Not yet tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_read_binarys(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_read_binaryd(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_read_binarysc(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_read_binarydc(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_read_binaryi(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_read_binaryl(iUnit,adata,nrec,ioprec,lchar,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_io
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
      INTEGER                          , INTENT(IN   ) :: iUnit,nrec,ioprec
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(  OUT) :: adata
#endif
      INTEGER                          , INTENT(  OUT) :: info
      LOGICAL                          , INTENT(IN   ) :: lchar
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0
      INTEGER                          :: ndata,iopt
      CHARACTER(LEN=ppm_char)          :: mesg
      INTEGER, DIMENSION(1)            :: ldl
#if __KIND == __INTEGER
      CHARACTER(LEN=SIZE(adata,1)), POINTER :: cbuf
      INTEGER                          :: ilen,i,ishift
#endif

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_read_binary',t0,info)
      ndata = SIZE(adata,1)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
      ENDIF

      !-------------------------------------------------------------------------
      !  Check for read-past-end-of-record exception
      !-------------------------------------------------------------------------
      IF (nrec .LT. ndata) THEN
          info = ppm_error_error
          WRITE(mesg,'(A,I11,A,I11)') 'Record lendth: ',nrec,   &
     &        '; data size: ',ndata
          CALL ppm_error(ppm_err_read_eor,'ppm_io_read_binary',mesg,   &
     &        __LINE__,info)
          GOTO 9999
      ENDIF

#if __KIND == __INTEGER
      !-------------------------------------------------------------------------
      !  Allocate character buffer if needed
      !-------------------------------------------------------------------------
      IF (lchar) THEN
          ALLOCATE(cbuf,STAT=info)
          IF (info .NE. 0) GOTO 800
      ENDIF
#endif
      
      !-------------------------------------------------------------------------
      !  Read actual data
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      IF (ioprec .EQ. ppm_param_io_double) THEN
          iopt = ppm_param_alloc_fit
          ldl(1) = ndata
          CALL ppm_alloc(dbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_io_read_binary',    &
     &            'double buffer DBUFFER',__LINE__,info)
              GOTO 9999
          ENDIF
          READ(iUnit,ERR=800,IOSTAT=info) dbuffer
          adata = REAL(dbuffer,ppm_kind_single)
          iopt  = ppm_param_dealloc
          CALL ppm_alloc(dbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_read_binary',    &
     &            'double buffer DBUFFER',__LINE__,info)
          ENDIF
      ELSE
          READ(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      IF (ioprec .EQ. ppm_param_io_single) THEN
          iopt = ppm_param_alloc_fit
          ldl(1) = ndata
          CALL ppm_alloc(sbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_io_read_binary',    &
     &            'single buffer SBUFFER',__LINE__,info)
              GOTO 9999
          ENDIF
          READ(iUnit,ERR=800,IOSTAT=info) sbuffer
          adata = REAL(sbuffer,ppm_kind_double)
          iopt  = ppm_param_dealloc
          CALL ppm_alloc(sbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_read_binary',    &
     &            'single buffer SBUFFER',__LINE__,info)
          ENDIF
      ELSE
          READ(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      IF (ioprec .EQ. ppm_param_io_double) THEN
          iopt = ppm_param_alloc_fit
          ldl(1) = ndata
          CALL ppm_alloc(dcbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_io_read_binary',    &
     &            'double complex buffer DCBUFFER',__LINE__,info)
              GOTO 9999
          ENDIF
          READ(iUnit,ERR=800,IOSTAT=info) dcbuffer
          adata = CMPLX(dcbuffer,KIND=ppm_kind_single)
          iopt  = ppm_param_dealloc
          CALL ppm_alloc(dcbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_read_binary',    &
     &            'double complex buffer DCBUFFER',__LINE__,info)
          ENDIF
      ELSE
          READ(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (ioprec .EQ. ppm_param_io_single) THEN
          iopt = ppm_param_alloc_fit
          ldl(1) = ndata
          CALL ppm_alloc(scbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_io_read_binary',    &
     &            'single complex buffer SCBUFFER',__LINE__,info)
              GOTO 9999
          ENDIF
          READ(iUnit,ERR=800,IOSTAT=info) scbuffer
          adata = CMPLX(scbuffer,KIND=ppm_kind_double)
          iopt  = ppm_param_dealloc
          CALL ppm_alloc(scbuffer,ldl,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_read_binary',    &
     &            'single complex buffer SCBUFFER',__LINE__,info)
          ENDIF
      ELSE
          READ(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __INTEGER
      IF (lchar) THEN
          READ(iUnit,ERR=800,IOSTAT=info) cbuf
      ELSE
          READ(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#else
      READ(iUnit,ERR=800,IOSTAT=info) adata
#endif
      IF (info .NE. 0) GOTO 800

#if __KIND == __INTEGER
      !-------------------------------------------------------------------------
      !  Convert to integer array and deallocate character buffer
      !-------------------------------------------------------------------------
      IF (lchar) THEN
          ilen = SIZE(adata,1)
          ishift = IACHAR('a')
          DO i=1,ilen
              adata(i) = IACHAR(cbuf(i:i)) - ishift
          ENDDO
          DEALLOCATE(cbuf,STAT=info)
      ENDIF
#endif

      GOTO 9999
      !-------------------------------------------------------------------------
      !  Error handler
      !-------------------------------------------------------------------------
  800 CONTINUE
      info = ppm_error_error
      CALL ppm_error(ppm_err_io,'ppm_io_read_binary',    &
     &    'Error while reading from file',__LINE__,info)
      GOTO 9999

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_read_binary',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_read_binarys
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_read_binaryd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_read_binarysc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_read_binarydc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_read_binaryi
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_read_binaryl
#endif
