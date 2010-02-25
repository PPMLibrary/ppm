      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_map_connect
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is an interface for the connection
      !                 handling, say distributing and pruning the connections
      !
      !  Input        : lda       (I) : leading dimension of cd(:,:)
      !                 id(:)     (I) : local to global particle number mapping
      !                 Npart     (I) : number of particles
      !                 contype   (I) : action type
      !                                   ppm_param_connect_distribute
      !                                   ppm_param_connect_send
      !                                   ppm_param_connect_prune
      !                 lsymm     (L) : using symmetry or not (only used for
      !                                 pruning)
      !
      !  Input/output : cd(:,:)   (I) : connection data
      !                 Ncon      (I) : number of connections
      !
      !  Output       : info      (I) : Return status. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_connect.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2004/10/01 16:33:35  ivos
      !  cosmetics.
      !
      !  Revision 1.4  2004/10/01 16:09:03  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.3  2004/07/26 11:46:54  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 08:54:23  ivos
      !  Changed names of renamed subroutines.
      !
      !  Revision 1.1  2004/07/26 08:24:43  ivos
      !  Check-in after module cleanup. All ppm_map_connect* were
      !  renamed to ppm_map_connect*. Also the modules.
      !
      !  Revision 1.8  2004/07/26 07:42:43  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/07/23 17:36:31  oingo
      !  Changed id from POINTER to INTENT(IN)
      !
      !  Revision 1.6  2004/07/20 08:58:37  oingo
      !  The id array is now scalar (i.e. dummy leading dimension removed)
      !
      !  Revision 1.5  2004/07/16 14:46:28  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.4  2004/07/16 14:07:59  ivos
      !  Added sequence and argument checks. New checks now allow multiple
      !  push-send-pop cycles per mapping.
      !
      !  Revision 1.3  2004/07/14 11:49:12  oingo
      !  Cosmetics
      !
      !  Revision 1.2  2004/06/10 16:20:00  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/05/11 13:26:51  oingo
      !  Initial release
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_connect(cd,lda,Ncon,id,Npart,contype,lsymm,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_map_connect_distrib
      USE ppm_module_map_connect_send
      USE ppm_module_map_connect_prune
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), POINTER       :: cd
      INTEGER                , INTENT(IN   ) :: lda
      INTEGER                , INTENT(INOUT) :: Ncon
      INTEGER, DIMENSION(:)  , INTENT(IN   ) :: id
      INTEGER                , INTENT(IN   ) :: Npart
      INTEGER                , INTENT(IN   ) :: contype
      LOGICAL                , INTENT(IN   ) :: lsymm
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_connect',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_connect',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Ncon .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_connect', &
     &            'Number of connections must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_connect',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_connect',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Distribute the connections to all processors
      !-------------------------------------------------------------------------
      IF (contype .EQ. ppm_param_connect_distribute) THEN
          CALL ppm_map_connect_distrib(cd,lda,Ncon,id,Npart,info)
          IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Distribute the connectins according to the current mapping
      !-------------------------------------------------------------------------
      ELSE IF (contype .EQ. ppm_param_connect_send) THEN
          CALL ppm_map_connect_send(cd,lda,Ncon,id,Npart,info)
          IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Prune the connections that are not needed
      !-------------------------------------------------------------------------
      ELSE IF (contype .EQ. ppm_param_connect_prune) THEN
          CALL ppm_map_connect_prune(cd,lda,Ncon,id,Npart,lsymm,info)
          IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Unknown connection mapping
      !-------------------------------------------------------------------------
      ELSE
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_map_connect',   &
     &        'Unknown connection mapping',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_connect',t0,info)
      RETURN

      END SUBROUTINE ppm_map_connect
