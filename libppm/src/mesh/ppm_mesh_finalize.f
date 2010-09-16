      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_finalize
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_finalize(info)
      !!! This routine terminates the mesh part of the ppm
      !!! library by deallocating all mesh memory structures.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
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
