      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_topo_dealloc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine deallocates all members of a 
      !                 topology object.
      !
      !  Input        : 
      !
      !  Input/output : topo  (ppm_t_topo) Topology structure to be
      !                                deallocated
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Remarks      : 
      !                
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_dealloc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:08  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_topo_dealloc(topo,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_typedef
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo)        , INTENT(INOUT) :: topo
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3) :: ldc
      INTEGER                :: iopt
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_dealloc',t0,info)
      ldc  = 1
      iopt = ppm_param_dealloc

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate all array members
      !-------------------------------------------------------------------------
      CALL ppm_alloc(topo%min_subs,ldc,iopt,info)
      CALL ppm_alloc(topo%max_subs,ldc,iopt,info)
      CALL ppm_alloc(topo%min_subd,ldc,iopt,info)
      CALL ppm_alloc(topo%max_subd,ldc,iopt,info)
      CALL ppm_alloc(topo%subs2proc,ldc,iopt,info)
      CALL ppm_alloc(topo%isublist,ldc,iopt,info)
      CALL ppm_alloc(topo%bcdef,ldc,iopt,info)
      CALL ppm_alloc(topo%sub_costs,ldc,iopt,info)
      CALL ppm_alloc(topo%sub_costd,ldc,iopt,info)
      CALL ppm_alloc(topo%nneighsubs,ldc,iopt,info)
      CALL ppm_alloc(topo%ineighsubs,ldc,iopt,info)
      CALL ppm_alloc(topo%subs_bc,ldc,iopt,info)
      CALL ppm_alloc(topo%ineighproc,ldc,iopt,info)
      CALL ppm_alloc(topo%icommseq,ldc,iopt,info)
      CALL ppm_alloc(topo%min_physs,ldc,iopt,info)
      CALL ppm_alloc(topo%min_physd,ldc,iopt,info)
      CALL ppm_alloc(topo%max_physs,ldc,iopt,info)
      CALL ppm_alloc(topo%max_physd,ldc,iopt,info)

      !-------------------------------------------------------------------------
      !  Deallocate the meshes
      !-------------------------------------------------------------------------
      CALL ppm_mesh_alloc_equi(topo%mesh,ldc,iopt,info)
      NULLIFY(topo%mesh)
      topo%max_meshid = 0

      !-------------------------------------------------------------------------
      !  Mark this topology as not defined
      !-------------------------------------------------------------------------
      topo%isdefined = .FALSE.

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_dealloc',t0,info)
      RETURN
      END SUBROUTINE ppm_topo_dealloc
