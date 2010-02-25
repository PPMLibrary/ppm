      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_connect_prune
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine prunes all connections that are not needed
      !                 on this processor according to the current topology
      !                 and whether symmetry will be used or not
      !
      !  Input        : lda       (I) : leading dimension of cd(:,:)
      !                 id(:)     (I) : local to global particle number mapping
      !                 Npart     (I) : number of particles
      !                 lsymm     (L) : symmetry or not
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
      !  $Log: ppm_map_connect_prune.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2004/11/11 15:24:33  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.4  2004/10/01 16:33:35  ivos
      !  cosmetics.
      !
      !  Revision 1.3  2004/10/01 16:09:03  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.2  2004/07/26 08:54:24  ivos
      !  Changed names of renamed subroutines.
      !
      !  Revision 1.1  2004/07/26 08:24:43  ivos
      !  Check-in after module cleanup. All ppm_map_part_connect* were
      !  renamed to ppm_map_connect*. Also the modules.
      !
      !  Revision 1.6  2004/07/26 07:42:44  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.5  2004/07/23 17:39:22  oingo
      !  Changed id from POINTER to INTENT(IN)
      !
      !  Revision 1.4  2004/07/20 08:59:18  oingo
      !  The id array is now scalar (i.e. dummy leading dimension removed)
      !
      !  Revision 1.3  2004/07/14 11:52:41  oingo
      !  Cosmetics
      !
      !  Revision 1.2  2004/06/10 16:20:01  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/05/11 13:28:01  oingo
      !  Initial release
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_connect_prune(cd,lda,Ncon,id,Npart,lsymm,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_invert_list
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
      LOGICAL                , INTENT(IN   ) :: lsymm
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: i,j,k,iopt,icon,part,prunemode
      REAL(ppm_kind_double) :: t0
      INTEGER               :: ncons_local

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_connect_prune',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Ncon .LT. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_argument,      &
     &            'ppm_map_connect_prune',     &
     &            'Number of connections must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we have a ring topology ...
      !-------------------------------------------------------------------------
      prunemode = 0
      IF (ppm_topoid .EQ. 0) THEN
          prunemode = 1
      !-------------------------------------------------------------------------
      !  ... or any other valid topology
      !-------------------------------------------------------------------------
      ELSE IF (ppm_topoid .GT. 0) THEN
          IF (lsymm) THEN
              ! if all particles are there
              prunemode = 2
          ELSE
              ! lowest id
              prunemode = 1
          ENDIF
      ELSE
          ! no topology defined
          prunemode = 0
      ENDIF

      IF (prunemode .NE. 0) THEN
          !---------------------------------------------------------------------
          !  Allocate
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = lda
          ldu(2) = Ncon
          CALL ppm_alloc(cd_local,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_connect_prune',  &
     &            'allocating connection list CD_LOCAL',__LINE__,info)
              GOTO 9999
          ENDIF

          cd_local(:,:) = cd(:,:)
          ncons_local = Ncon

          !---------------------------------------------------------------------
          !  Get the inverse of the local to global particle id mapping
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = Npart
          CALL ppm_alloc(id_temp,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_connect_prune',  &
     &            'allocating connection array ID_TEMP',__LINE__,info)
              GOTO 9999
          ENDIF

          id_temp(1:Npart) = id(1:Npart)

          CALL ppm_util_invert_list(id_temp,id_inv,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_connect_prune',  &
     &            'inverting id list',__LINE__,info)
              GOTO 9999
          ENDIF

          iopt = ppm_param_dealloc
          CALL ppm_alloc(id_temp,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_connect_prune',  &
     &            'deallocating connection array ID_TEMP',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  prunemode 1: remove all connections where we do not have the particle
      !               with the lowest ID of the connection
      !-------------------------------------------------------------------------
      IF (prunemode .EQ. 1) THEN
          !---------------------------------------------------------------------
          !  Loop over all connections
          !---------------------------------------------------------------------
          icon = 0
          DO j = 1,ncons_local
             !------------------------------------------------------------------
             !  Get the particle with the lowest ID in this connection
             !------------------------------------------------------------------
             part = MINVAL(cd_local(:,j))

             !------------------------------------------------------------------
             !  Check if we have this particle
             !------------------------------------------------------------------
             IF ((part .GE. LBOUND(id_inv,1)) .AND.  &
     &           (part .LE. UBOUND(id_inv,1))) THEN
                 IF (id_inv(part) .EQ. -HUGE(part)) THEN
                     cd_local(1,j) = -1
                     icon = icon + 1
                 ENDIF
             ELSE
                 cd_local(1,j) = -1
                 icon = icon + 1
             ENDIF
          ENDDO
      !-------------------------------------------------------------------------
      !  prunemode 2: remove all connections we do not have all particles
      !               for
      !-------------------------------------------------------------------------
      ELSE IF (prunemode .EQ. 2) THEN
          icon = 0
          DO j = 1,ncons_local
             k = 0
             !------------------------------------------------------------------
             !  Loop over all particles in the connection and count how
             !  many are missing.
             !------------------------------------------------------------------
             DO i = 1,lda
                IF ((cd_local(i,j) .GE. LBOUND(id_inv,1)) .AND.  &
     &              (cd_local(i,j) .LE. UBOUND(id_inv,1))) THEN
                    IF (id_inv(cd_local(i,j)) .EQ. -HUGE(cd_local(i,j))) THEN
                       k = k + 1
                    ENDIF
                ELSE
                    k = k + 1
                ENDIF
             ENDDO

             !------------------------------------------------------------------
             !  If at least one particle is missing we prune this connection
             !------------------------------------------------------------------
             IF (k .GT. 0) THEN
                 cd_local(1,j) = -1
                 icon = icon + 1
             ENDIF
          ENDDO
      ENDIF

      IF (prunemode .NE. 0) THEN
          !---------------------------------------------------------------------
          !  Prepare the final connection array
          !---------------------------------------------------------------------
          ldu(1) = lda
          ldu(2) = Ncon - icon
          iopt = ppm_param_alloc_fit
          CALL ppm_alloc(cd,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_connect_prune',  &
     &            'allocating connection array CD',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Go through all connections and save all that are not marked
          !---------------------------------------------------------------------
          Ncon = 0
          DO j = 1,ncons_local
             IF (cd_local(1,j) .NE. -1) THEN
                 Ncon = Ncon + 1
                 cd(:,Ncon) = cd_local(:,j)
             ENDIF
          ENDDO

          iopt = ppm_param_dealloc
          CALL ppm_alloc(cd_local,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_connect_prune',  &
     &            'deallocating connection list CD_LOCAL',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_connect_prune',t0,info)
      RETURN

      END SUBROUTINE ppm_map_connect_prune
