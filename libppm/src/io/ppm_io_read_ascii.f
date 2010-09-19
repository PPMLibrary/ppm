      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_read_ascii
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_readascii_s(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_readascii_d(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE io_readascii_sc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE io_readascii_dc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_readascii_i(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE io_readascii_l(iUnit,adata,iofmt,lchar,info)
#endif
      !!! Provides the atomic reading functionality for ASCII 1D arrays.
      !!!
      !!! NOTE: Never used by the user program directly.

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
      !!! Fortran I/O unit to be used
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(  OUT) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(  OUT) :: adata
#endif
      !!! Data vector to be read. 1d array of either single, double,
      !!! single complex, double complex, integer or logical elements.
      !!! Needs to be allocated to proper size before calling this
      !!! routine. Its size determines the number of elements read from the
      !!! file.
      INTEGER                          , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      CHARACTER(LEN=*)                 , INTENT(IN   ) :: iofmt
      !!! Format string for I/O. List-directed I/O is used if
      !!! length zero string is given.
      LOGICAL                          , INTENT(IN   ) :: lchar
      !!! `.TRUE.` if a CHARACTER string is to be read and has to be
      !!! converted to an INTEGER array before returning it.
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
      END SUBROUTINE io_readascii_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_readascii_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE io_readascii_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE io_readascii_dc
#elif __KIND == __INTEGER
      END SUBROUTINE io_readascii_i
#elif __KIND == __LOGICAL
      END SUBROUTINE io_readascii_l
#endif
