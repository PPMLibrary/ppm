      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_write_ascii
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
      SUBROUTINE io_writeascii_s(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_writeascii_d(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE io_writeascii_sc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE io_writeascii_dc(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_writeascii_i(iUnit,adata,iofmt,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE io_writeascii_l(iUnit,adata,iofmt,lchar,info)
#endif
      !!! Provides the atomic writing functionality for ascii 1D arrays.
      !!!
      !!! NOTE: It is never used by the user program directly.

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
      REAL(MK)   , DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(IN   ) :: adata
#endif
      !!! Data vector to be read. 1d array of either single, double,
      !!! single complex, double complex, integer or logical elements. All
      !!! elements of it will be written.
      INTEGER                          , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      CHARACTER(LEN=*)                 , INTENT(IN   ) :: iofmt
      !!! Format string for I/O. List-directed I/O is used if
      !!! length zero string is given.
      LOGICAL                          , INTENT(IN   ) :: lchar
      !!! .TRUE. if the integer array actually is a CHARACTER string
      !!! and should be converted back before output.
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
      END SUBROUTINE io_writeascii_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_writeascii_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE io_writeascii_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE io_writeascii_dc
#elif __KIND == __INTEGER
      END SUBROUTINE io_writeascii_i
#elif __KIND == __LOGICAL
      END SUBROUTINE io_writeascii_l
#endif
