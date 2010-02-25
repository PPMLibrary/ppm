      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_define
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine defines a new mesh on an existing
      !                 topology. The routine checks that the subdomains of
      !                 the topology are compatible with the specified mesh
      !                 (i.e. are integer multiples of the mesh spacing in
      !                 extent). If not, an error is returned.
      !
      !  Input        : topo_id      (I) topology ID for which to create
      !                                  mesh (user numbering)
      !                 Nm(:)        (I) number of mesh POINTS in each
      !                                  dimension. Subs must be compatible
      !                                  with this mesh, otherwise an error
      !                                  occurs.
      !                 min_phys(:)  (F) min. extent of the computational
      !                                  domain
      !                 max_phys(:)  (F) max. extent of the computational
      !                                  domain
      !
      !  Input/output : mesh_id      (I) mesh ID of the new meshi (user
      !                                  numbering). 
      !                                  If .LE. 0 on input, the 
      !                                  routine will create an automatic 
      !                                  one and return it here.
      !
      !  Output       : istart(:,:)  (I) start indices of all subs meshes
      !                                  in global mesh
      !                 ndata(:,:)   (I) number of mesh points in each
      !                                  direction on each sub mesh
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mesh_define.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2004/10/01 16:09:09  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.7  2004/08/31 13:29:58  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.6  2004/07/26 13:49:17  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.5  2004/07/26 11:48:09  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.4  2004/07/26 07:42:48  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.3  2004/07/16 14:46:28  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.2  2004/04/22 10:35:55  ivos
      !  bugfix: use of SIZE(ppm_internal_topid) caused problems in argument
      !  checking when first user-topoid was > 1. Resolved by use of UBOUND.
      !
      !  Revision 1.1  2004/03/06 17:24:03  ivos
      !  Initial implementation. Tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mesh_define_s(topo_id,Nm,min_phys,max_phys,mesh_id, &
     &                             istart,ndata,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mesh_define_d(topo_id,Nm,min_phys,max_phys,mesh_id, &
     &                             istart,ndata,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_mesh_on_subs
      USE ppm_module_mesh_store
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topo_id
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys,max_phys
      INTEGER                 , INTENT(INOUT) :: mesh_id
      INTEGER                 , INTENT(  OUT) :: info
      INTEGER , DIMENSION(:,:), POINTER       :: istart,ndata
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: t0
      INTEGER                :: topoid
      LOGICAL                :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_define',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_mesh_define',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_user,topo_id,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_define',  &
     &            'topo_id is invalid!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Translate user topo_id to internal topoid
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)

      !-------------------------------------------------------------------------
      !  Create new mesh
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,ppm_min_subs(:,:,topoid),  &
     &    ppm_max_subs(:,:,topoid),ppm_nsubs(topoid),istart,ndata,info)
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,ppm_min_subd(:,:,topoid),  &
     &    ppm_max_subd(:,:,topoid),ppm_nsubs(topoid),istart,ndata,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_error
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new mesh
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(mesh_id,topoid,ppm_nsubs(topoid),ndata,istart,Nm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_mesh_define',   &
     &        'Storing new mesh failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_define',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mesh_define_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mesh_define_d
#endif
