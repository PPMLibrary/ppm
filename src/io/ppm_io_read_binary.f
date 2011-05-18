      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_io_read_binary
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
      SUBROUTINE io_readbin_s(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_readbin_d(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE io_readbin_sc(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE io_readbin_dc(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_readbin_i(iUnit,adata,nrec,ioprec,lchar,info)
#elif __KIND == __LOGICAL
      SUBROUTINE io_readbin_l(iUnit,adata,nrec,ioprec,lchar,info)
#endif
      !!! Provides the atomic reading functionality for binary 1D arrays.
      !!!
      !!! NOTE: It is never used by the user program directly.

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
      INTEGER                          , INTENT(IN   ) :: iUnit
      !!! Fortran I/O unit to be used
      INTEGER                          , INTENT(IN   ) :: nrec
      !!! Length of record to be read
      INTEGER                          , INTENT(IN   ) :: ioprec
      !!! Precision for I/O. One of:
      !!!
      !!! * ppm_param_io_single
      !!! * ppm_param_io_double
      !!! 
      !!! Integer and Logical cases are not influenced by this. No
      !!! conversions are done if any value other than the above is given.
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
      LOGICAL                          , INTENT(IN   ) :: lchar
      !!! `.TRUE.` if a CHARACTER string is to be read and has to be
      !!! converted to an INTEGER array before returning it.
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
      END SUBROUTINE io_readbin_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_readbin_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE io_readbin_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE io_readbin_dc
#elif __KIND == __INTEGER
      END SUBROUTINE io_readbin_i
#elif __KIND == __LOGICAL
      END SUBROUTINE io_readbin_l
#endif
