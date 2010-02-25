      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_io_write_binary
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine provides the atomic writing
      !                 functionality for binary 1d arrays.
      !                 It is never used by the user program directly.
      !
      !  Input        : iUnit          (I) Fortran I/O unit to be used
      !                 ioprec         (I) Precision for I/O.
      !                                    One of:
      !                                       ppm_param_io_single
      !                                       ppm_param_io_double
      !                                    Integer and Logical cases are
      !                                    not influenced by this. No
      !                                    conversions are done if any
      !                                    value other than the above is
      !                                    given.
      !                 adata          (O) Data vector to be written. 1d
      !                                    array of either single, double,
      !                                    single complex, double complex,
      !                                    integer or logical elements. All
      !                                    elements of it will be written.
      !                 lchar          (L) .TRUE. if the integer array
      !                                    actually is a CHARACTER string
      !                                    and should be converted back
      !                                    before output.
      !                                    
      !  Input/output : 
      !
      !  Output       : info           (I) return code. 0 on success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_write_binary.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2006/02/03 09:35:19  ivos
      !  Fixed bug 00004: cbuf was only allocated in the __INTEGER case, but
      !  used in an else-branch. Added separate ifdefs for the __INTEGER
      !  case to make sure cbuf is ONLY used there.
      !
      !  Revision 1.6  2004/10/01 16:09:02  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/07/26 07:45:28  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/06/11 10:21:58  ivos
      !  Changed reclen from I8B=SELECTED_INT_KIND(18) to regular INTEGER
      !  since I8B caused problems on the NEC.
      !
      !  Revision 1.3  2004/06/04 08:43:16  ivos
      !  Corrected typo in error message.
      !
      !  Revision 1.2  2004/05/27 14:52:11  ivos
      !  Bugfix: fixed misplaced GOTO
      !
      !  Revision 1.1  2004/05/27 10:43:41  ivos
      !  Initial implementation. Not yet tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_io_write_binarys(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_write_binaryd(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_write_binarysc(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_write_binarydc(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_write_binaryi(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_write_binaryl(iUnit,adata,ioprec,lchar,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      INTEGER                          , INTENT(IN   ) :: iUnit,ioprec
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(IN   ) :: adata
#endif
      INTEGER                          , INTENT(  OUT) :: info
      LOGICAL                          , INTENT(IN   ) :: lchar
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0
      INTEGER                          :: reclen
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
      CALL substart('ppm_io_write_binary',t0,info)
      reclen = SIZE(adata,1)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
      ENDIF

#if __KIND == __INTEGER
      !-------------------------------------------------------------------------
      !  Convert to character string if needed
      !-------------------------------------------------------------------------
      IF (lchar) THEN
          ALLOCATE(cbuf,STAT=info)
          IF (info .NE. 0) GOTO 800
          ilen = SIZE(adata,1)
          ishift = IACHAR('a')
          DO i=1,ilen
              cbuf(i:i) = CHAR(adata(i)+ishift)
          ENDDO
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Write record length to file
      !-------------------------------------------------------------------------
      WRITE(iUnit,ERR=800,IOSTAT=info) reclen
      IF (info .NE. 0) GOTO 800

      !-------------------------------------------------------------------------
      !  Write actual data
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION 
      IF (ioprec .EQ. ppm_param_io_double) THEN
          WRITE(iUnit,ERR=800,IOSTAT=info) REAL(adata,ppm_kind_double)
      ELSE
          WRITE(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __DOUBLE_PRECISION 
      IF (ioprec .EQ. ppm_param_io_single) THEN
          WRITE(iUnit,ERR=800,IOSTAT=info) REAL(adata,ppm_kind_single)
      ELSE
          WRITE(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __SINGLE_PRECISION_COMPLEX 
      IF (ioprec .EQ. ppm_param_io_double) THEN
          WRITE(iUnit,ERR=800,IOSTAT=info) CMPLX(adata,KIND=ppm_kind_double)
      ELSE
          WRITE(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __DOUBLE_PRECISION_COMPLEX 
      IF (ioprec .EQ. ppm_param_io_single) THEN
          WRITE(iUnit,ERR=800,IOSTAT=info) CMPLX(adata,KIND=ppm_kind_single)
      ELSE
          WRITE(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#elif __KIND == __INTEGER
      IF (lchar) THEN
          WRITE(iUnit,ERR=800,IOSTAT=info) cbuf
      ELSE
          WRITE(iUnit,ERR=800,IOSTAT=info) adata
      ENDIF
#else
      WRITE(iUnit,ERR=800,IOSTAT=info) adata
#endif
      IF (info .NE. 0) GOTO 800

#if __KIND == __INTEGER
      !-------------------------------------------------------------------------
      !  Deallocate character buffer if needed
      !-------------------------------------------------------------------------
      IF (lchar) THEN
          DEALLOCATE(cbuf,STAT=info)
      ENDIF
#endif

      GOTO 9999
      !-------------------------------------------------------------------------
      !  Error handler
      !-------------------------------------------------------------------------
  800 CONTINUE
      info = ppm_error_error
      CALL ppm_error(ppm_err_io,'ppm_io_write_binary',    &
     &    'Error while writing to file',__LINE__,info)
      GOTO 9999

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_write_binary',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_write_binarys
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_write_binaryd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_write_binarysc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_write_binarydc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_write_binaryi
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_write_binaryl
#endif
