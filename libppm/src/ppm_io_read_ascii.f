      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_read_ascii
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine provides the atomic reading
      !                 functionality for ascii 1d arrays.
      !                 It is never used by the user program directly.
      !
      !  Input        : iUnit          (I) Fortran I/O unit to be used
      !                 iofmt          (C) Format string for I/O.
      !                                    List-directed I/O is used if
      !                                    length zero string is given.
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
      !  $Log: ppm_io_read_ascii.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2006/02/03 09:35:20  ivos
      !  Fixed bug 00004: cbuf was only allocated in the __INTEGER case, but
      !  used in an else-branch. Added separate ifdefs for the __INTEGER
      !  case to make sure cbuf is ONLY used there.
      !
      !  Revision 1.4  2004/10/01 16:09:01  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
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
      SUBROUTINE ppm_io_read_asciis(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_io_read_asciid(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_read_asciisc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_io_read_asciidc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_io_read_asciii(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_io_read_asciil(iUnit,adata,iofmt,lchar,info)
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
      REAL(MK)   , DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(  OUT) :: adata
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
      CALL substart('ppm_io_read_ascii',t0,info)

      IF (LEN_TRIM(iofmt) .EQ. 0) THEN
          lfmt = .FALSE.
      ELSE
          lfmt = .TRUE.
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
      IF (lfmt) THEN
#if __KIND == __INTEGER
          IF (lchar) THEN
              READ(iUnit,iofmt,ERR=800,IOSTAT=info) cbuf
          ELSE
              READ(iUnit,iofmt,ERR=800,IOSTAT=info) adata
          ENDIF
#else
          READ(iUnit,iofmt,ERR=800,IOSTAT=info) adata
#endif
      ELSE
#if __KIND == __INTEGER
          IF (lchar) THEN
              READ(iUnit,*,ERR=800,IOSTAT=info) cbuf
          ELSE
              READ(iUnit,*,ERR=800,IOSTAT=info) adata
          ENDIF
#else
          READ(iUnit,*,ERR=800,IOSTAT=info) adata
#endif
      ENDIF
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
      CALL ppm_error(ppm_err_io,'ppm_io_read_ascii',    &
     &    'Error while reading from file',__LINE__,info)
      GOTO 9999

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_read_ascii',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_io_read_asciis
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_io_read_asciid
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_read_asciisc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_io_read_asciidc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_io_read_asciii
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_io_read_asciil
#endif
