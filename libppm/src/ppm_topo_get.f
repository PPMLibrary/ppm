      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_get
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Reads the information of an internally stored
      !                 topology and returns it to the user.
      !
      !  Input        : topo_id    (I) ID of the topology to be read.
      !
      !  Output       : topo   (ppm_t_topo) A copy of the topology object
      !                                as stored internally.
      !                 info       (I) return status. 
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_get.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:08  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_topo_get(topo_id,topo,info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_check_topoid
      USE ppm_module_typedef
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(  OUT) :: info
      INTEGER                 , INTENT(IN   ) :: topo_id
      TYPE(ppm_t_topo)        , INTENT(  OUT) :: topo
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)     :: t0
      LOGICAL                   :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_get',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_get',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_user,topo_id,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_get',  &
     &            'Topology ID is invalid!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Read the topology
      !-------------------------------------------------------------------------
      topo%ID = topo_id

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_get',t0,info)
      RETURN
      END SUBROUTINE ppm_topo_get
