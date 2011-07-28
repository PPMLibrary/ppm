      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_remap
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps the particles onto the given topology.
      !
      !  Input        : xp(:,:)    (F) the position of the particles
      !                 Npart      (I) the number of particles 
      !                                (on processor)
      !                 topoid     (I) topology identifier (internal
      !                                numbering) of destination topology.
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_remap.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2004/08/31 12:48:10  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.7  2004/07/26 07:42:46  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.6  2004/06/10 16:20:02  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.5  2004/01/23 17:24:17  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.4  2004/01/23 11:31:23  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.3  2004/01/06 13:53:13  ivos
      !  Changed INTENT of topoid to IN, removed 2 non-used local 
      !  variables, removed non-used #ifdef at the end and switched use of 
      !  print to pwrite.
      !
      !  Revision 1.2  2003/12/17 15:47:45  hiebers
      !  new implementation
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_remap_s(xp,Npart,topoid,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_remap_d(xp,Npart,topoid,info)
#endif 
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
      USE ppm_module_topo_check
      USE ppm_module_check_topoid
      USE ppm_module_map_part_global
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:)  , INTENT(IN   ) :: xp
      INTEGER                   , INTENT(IN   ) :: Npart
      INTEGER                   , INTENT(IN   ) :: topoid
      INTEGER                   , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      LOGICAL             :: topo_ok
      REAL(MK)            :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_remap',t0,info)
      
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Npart .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_remap',  &
     &            'Npart must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,topo_ok,info)
          IF (.NOT. topo_ok) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_remap',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if current topology is ok
      !-------------------------------------------------------------------------
      CALL ppm_topo_check(xp,Npart,topo_ok,info)

      !-------------------------------------------------------------------------
      !  if topology is ok remap particles otherwise abort remapping
      !-------------------------------------------------------------------------
      IF(topo_ok) THEN
         CALL ppm_map_part_global(xp,Npart,topoid,info)
      ELSE
         info = ppm_error_error
         CALL ppm_error(ppm_err_topo_missm,'ppm_map_part_remap',  &
     &     'Particles are not on current topology. Mapping discarded.',  &
     &     __LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_remap',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_remap_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_remap_d
#endif
