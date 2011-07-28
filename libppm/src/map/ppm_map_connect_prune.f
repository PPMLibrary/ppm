      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_connect_prune
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_connect_prune(topoid,cd,lda,Ncon,id,Npart,lsymm,info)
      !!! This routine prunes all connections that are not needed
      !!! on this processor according to the current topology
      !!! and whether symmetry will be used or not.

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
      INTEGER                , INTENT(IN   ) :: topoid
      !!! ID of current topology
      INTEGER, DIMENSION(:,:), POINTER       :: cd
      !!! Connection data
      INTEGER                , INTENT(IN   ) :: lda
      !!! Leading dimension of cd(:,:)
      INTEGER                , INTENT(INOUT) :: Ncon
      !!! Number of connections
      INTEGER, DIMENSION(:)  , INTENT(IN   ) :: id
      !!! Local to global particle number mapping
      INTEGER                , INTENT(IN   ) :: Npart
      !!! Number of particles
      LOGICAL                , INTENT(IN   ) :: lsymm
      !!! Is symmetric?
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: i,j,k,iopt,icon,part,prunemode
      REAL(ppm_kind_double) :: t0
      LOGICAL               :: valid
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
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we have a ring topology ...
      !-------------------------------------------------------------------------
      prunemode = 0
      IF (topoid .EQ. ppm_param_topo_undefined) THEN
          prunemode = 1
      !-------------------------------------------------------------------------
      !  ... or any other valid topology
      !-------------------------------------------------------------------------
      ELSE IF (topoid .NE. ppm_param_topo_undefined) THEN
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
      CONTAINS
      SUBROUTINE check
        IF (Ncon .LT. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_argument,      &
     &            'ppm_map_connect_prune',     &
     &            'Number of connections must be >=0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .NE. ppm_param_topo_undefined) THEN
            CALL ppm_check_topoid(topoid,valid,info)
            IF (.NOT. valid) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_map_connect_prune',  &
     &               'topoid out of range',__LINE__,info)
                GOTO 8888
            ENDIF
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_connect_prune
