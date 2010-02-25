      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_write_ascii
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine provides the atomic writing
      !                 functionality for ascii 1d arrays.
      !                 It is never used by the user program directly.
      !
      !  Input        : iUnit          (I) Fortran I/O unit to be used
      !                 iofmt          (C) Format string for I/O.
      !                                    List-directed I/O is used if
      !                                    length zero string is given.
      !                 adata          (O) Data vector to be read. 1d
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
      !  $Log: ppm_io_write_ascii.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/02/03 09:35:19  ivos
      !  Fixed bug 00004: cbuf was only allocated in the __INTEGER case, but
      !  used in an else-branch. Added separate ifdefs for the __INTEGER
      !  case to make sure cbuf is ONLY used there.
      !
      !  Revision 1.5  2004/10/01 16:09:02  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/26 07:45:28  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.3  2004/06/04 08:43:16  ivos
      !  Corrected typo in error message.
      !
      !  Revision 1.2  2004/05/27 14:52:10  ivos
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
      SUBROUTINE ppm_io_write_asciis(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_write_asciid(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_write_asciisc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_write_asciidc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_write_asciii(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_write_asciil(iUnit,adata,iofmt,lchar,info)
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
      INTEGER                          , INTENT(IN   ) :: iUnit
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
      CHARACTER(LEN=*)                 , INTENT(IN   ) :: iofmt
      LOGICAL                          , INTENT(IN   ) :: lchar
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0
      LOGICAL                          :: lfmt
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
      CALL substart('ppm_io_write_ascii',t0,info)

      IF (LEN_TRIM(iofmt) .EQ. 0) THEN
          lfmt = .FALSE.
      ELSE
          lfmt = .TRUE.
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
      !  Write actual data
      !-------------------------------------------------------------------------
      IF (lfmt) THEN
#if __KIND == __INTEGER
          IF (lchar) THEN
              WRITE(iUnit,iofmt,ERR=800,IOSTAT=info) TRIM(cbuf)
          ELSE
              WRITE(iUnit,iofmt,ERR=800,IOSTAT=info) adata
          ENDIF
#else
          WRITE(iUnit,iofmt,ERR=800,IOSTAT=info) adata
#endif
      ELSE
#if __KIND == __INTEGER
          IF (lchar) THEN
              WRITE(iUnit,*,ERR=800,IOSTAT=info) TRIM(cbuf)
          ELSE
              WRITE(iUnit,*,ERR=800,IOSTAT=info) adata
          ENDIF
#else
          WRITE(iUnit,*,ERR=800,IOSTAT=info) adata
#endif
      ENDIF
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
      CALL ppm_error(ppm_err_io,'ppm_io_write_ascii',    &
     &    'Error while writing to file',__LINE__,info)
      GOTO 9999

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_write_ascii',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_write_asciis
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_write_asciid
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_write_asciisc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_write_asciidc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_write_asciii
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_write_asciil
#endif
