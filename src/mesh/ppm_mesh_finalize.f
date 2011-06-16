      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_finalize
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

      SUBROUTINE ppm_mesh_finalize(info)
      !!! This routine terminates the mesh part of the ppm
      !!! library by deallocating all mesh memory structures.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mesh_alloc
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(OUT)    :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1)   :: lda 
      INTEGER                 :: iopt
      INTEGER                 :: istat
      REAL(ppm_kind_double)   :: t0
      CHARACTER(LEN=ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_finalize',t0,info)
      lda(1) = 0

      !-------------------------------------------------------------------------
      !  Deallocate mesh structures 
      !-------------------------------------------------------------------------
      istat = 0
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ppm_mesh_isendfromsub,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_isendblkstart,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_isendblksize,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_irecvtosub,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_irecvblkstart,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_irecvblksize,lda,iopt,info)
      istat = istat + info

      IF (istat .NE. 0) THEN
          WRITE(mesg,'(A,I3,A)') 'for ',istat,   &
     &        ' mesh arrays. Possible memory leak.'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_finalize',mesg,__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_finalize',t0,info)

      END SUBROUTINE ppm_mesh_finalize
