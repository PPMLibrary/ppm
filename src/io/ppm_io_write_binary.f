      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_io_write_binary
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
      SUBROUTINE io_writebin_s(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_writebin_d(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE io_writebin_sc(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE io_writebin_dc(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_writebin_i(iUnit,adata,ioprec,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE io_writebin_l(iUnit,adata,ioprec,lchar,info)
#endif
      !!! Provides the atomic writing functionality for binary 1D arrays.
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
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                          , INTENT(IN   ) :: iUnit
      !!! Fortran I/O unit to be used
      INTEGER                          , INTENT(IN   ) :: ioprec
      !!! Precision for I/O. One of:
      !!!
      !!! * ppm_param_io_single
      !!! * ppm_param_io_double
      !!!
      !!! Integer and Logical cases are not influenced by this. No
      !!! conversions are done if any value other than the above is given.
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)        , INTENT(IN   ) :: adata
#elif __KIND == __LOGICAL
      LOGICAL    , DIMENSION(:)        , INTENT(IN   ) :: adata
#endif
      !!! Data vector to be written. 1d array of either single, double,
      !!! single complex, double complex, integer or logical elements. All
      !!! elements of it will be written.
      INTEGER                          , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      LOGICAL                          , INTENT(IN   ) :: lchar
      !!! .TRUE. if the integer array actually is a CHARACTER string
      !!! and should be converted back before output.
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
      END SUBROUTINE io_writebin_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_writebin_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE io_writebin_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE io_writebin_dc
#elif __KIND == __INTEGER
      END SUBROUTINE io_writebin_i
#elif __KIND == __LOGICAL
      END SUBROUTINE io_writebin_l
#endif
