      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine terminates the mesh part of the ppm
      !                 library by deallocating all mesh memory structures.
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : info    (I) 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mesh_finalize.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/09/04 18:34:52  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.5  2004/10/01 16:33:38  ivos
      !  cosmetics.
      !
      !  Revision 1.4  2004/10/01 16:09:10  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.3  2004/07/26 11:48:09  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:42:48  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.1  2004/07/16 14:11:17  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_finalize(info)
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
      CALL ppm_mesh_alloc(ppm_meshid,lda,iopt,info)
      istat = istat + info
      CALL ppm_mesh_alloc(ppm_cart_mesh,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_max_meshid,lda,iopt,info)
      istat = istat + info
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
      CALL ppm_alloc(ppm_mesh_ghost_fromsub,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_tosub,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_blkstart,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_blksize,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_blk,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_nsend,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_nrecv,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_recvtosub,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_recvblkstart,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_recvblksize,lda,iopt,info)
      istat = istat + info
      CALL ppm_alloc(ppm_mesh_ghost_recvblk,lda,iopt,info)
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
