      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_map_field
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps a field onto a topology.
      !
      !  Input        : [lda]       (I) the leading dimension of the data
      !                                 (only important for push and pop)
      !                                 in the case of vector fields. For
      !                                 scalar fields, omit this argument.
      !                 to_topo     (I) user topology ID of destination. If 
      !                                 is <=0, the current topology is used.
      !                 from_mesh   (I) user mesh ID of source. If <=0,
      !                                 the first mesh of the current
      !                                 topology is used
      !                 to_mesh     (I) user mesh ID of destination. If <=0,
      !                                 the first mesh of the to_topo is used
      !                 ghostsize(:)(I) number of ghost mesh points in each
      !                                 dimension.
      !                 maptype     (I) Mapping action requested. One of:
      !                                     ppm_param_map_global
      !                                     ppm_param_map_push
      !                                     ppm_param_map_send
      !                                     ppm_param_map_pop
      !                                     ppm_param_map_cancel
      !
      !  Input/output : fv(:,...)   (O) field values (at mesh node
      !                                 positions) to be pushed/popped.
      !                                 Overloaded data types: single,
      !                                 double, integer, logical, single
      !                                 complex, double complex. Can be
      !                                 either a 3d (2d scalar field), 4d
      !                                 (2d vector field or 3d scalar
      !                                 field) or 5d (3d vector field)
      !                                 array. 1st index for vector fields
      !                                 is lda. For scalar fields it
      !                                 directly starts with the 2/3 mesh
      !                                 (i,j[,k]). Last index is subdomain ID
      !                                 (all local subs on the proc.).
      !                 mask(:,...) (L) Logical mask. Only the mesh nodes
      !                                 for which this is .TRUE. will be
      !                                 mapped. OPTIONAL. If not given, all
      !                                 points are mapped. After the send,
      !                                 the new mask is returned here.
      !                                 First 3 (2 for 2d meshes) indices
      !                                 are mesh (i,j[,k]), last one is
      !                                 subid for all subs on the local
      !                                 processor.
      !
      !  Output       : info        (I) return status. 0 upon success.
      !
      !  Remarks      : The current topology for fields is not given by
      !                 ppm_topoid, but ppm_field_topoid.
      !
      !                 The mesh ghost mapping is in ppm_map_field_ghost.f
      !
      !                 The argument ghostsize is only used by
      !                 map_field_pop to (re)allocate the fields to the
      !                 proper size. 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.23  2005/03/10 01:40:30  ivos
      !  Empty buffer is now only reported for ppm_debug.GT.0 and is no
      !  longer logged to prevent huge log files.
      !
      !  Revision 1.22  2004/08/31 12:14:54  ivos
      !  Changed to use ppm_check_topoid and ppm_check_meshid to check validity
      !  of user-numbered IDs.
      !
      !  Revision 1.21  2004/08/31 09:06:10  ivos
      !  Cleaned bugfix and removed PRINTs (see feedback on bug 00021).
      !
      !  Revision 1.20  2004/08/18 14:59:32  michaebe
      !  cosmetics
      !
      !  Revision 1.19  2004/08/18 14:53:06  michaebe
      !  validity check for mesh_id (see bug #0000021)
      !
      !  Revision 1.18  2004/07/26 11:46:54  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.17  2004/07/26 07:42:40  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.16  2004/07/22 10:50:51  hiebers
      !  Fixed Bug in the arguments of ppm_map_field_pop call
      !
      !  Revision 1.15  2004/07/19 11:01:50  ivos
      !  Overloaded field push and pop operations for scalar fields and
      !  added all changes to the module and the interface routines.
      !
      !  Revision 1.14  2004/07/16 14:46:27  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.13  2004/07/16 14:07:58  ivos
      !  Added sequence and argument checks. New checks now allow multiple
      !  push-send-pop cycles per mapping.
      !
      !  Revision 1.12  2004/04/22 10:35:55  ivos
      !  bugfix: use of SIZE(ppm_internal_topid) caused problems in argument
      !  checking when first user-topoid was > 1. Resolved by use of UBOUND.
      !
      !  Revision 1.11  2004/04/14 15:08:50  ivos
      !  Added OPTIONAL argument mask to allow selective mapping of only a
      !  subset of mesh points. Different masks can be used for different
      !  fields on the same mesh. Not tested yet.
      !
      !  Revision 1.10  2004/04/07 09:21:23  ivos
      !  Added ghostsize to the argument list since map_field_pop now needs 
      !  this.
      !
      !  Revision 1.9  2004/04/05 12:00:47  ivos
      !  Added security checks on ppm_map_type and changed CALL to 
      !  ppm_map_field_pop to new argument list.
      !
      !  Revision 1.8  2004/04/01 14:11:28  ivos
      !  Removed all the ghost stuff. Is now in ppm_map_field_ghost.f.
      !
      !  Revision 1.7  2004/02/24 12:16:04  ivos
      !  map_field_pop now grows the field data array if needed instead of just
      !  checking its size and requiring the user to allocate it. This made all
      !  arguments change from INTENT(INOUT) to POINTER. Reason for the change:
      !  symmetry with the particle pop.
      !
      !  Revision 1.6  2004/02/24 08:46:54  ivos
      !  Added overloading for single complex and double complex data types.
      !
      !  Revision 1.5  2004/02/20 16:26:01  ivos
      !  bugfix: ppm_nrecvbuffer was wrong in pop. and target_meshid was reset
      !  by send even though pop still needed it.
      !
      !  Revision 1.4  2004/02/20 15:45:51  ivos
      !  Added check if current (source) topology is defined at all.
      !
      !  Revision 1.3  2004/02/19 15:45:19  walther
      !  bug fix: added ppm_target_topoid, ppm_source_meshid, and 
      !  ppm_target_meshid.
      !
      !  Revision 1.2  2004/02/18 15:04:18  walther
      !  Renamed the ppm_param_map_ghost to ppm_param_map_ghost_put
      !
      !  Revision 1.1  2004/02/11 14:34:11  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      ! DIM is the data dimension (scalar or vector field)
      ! MESH_DIM is the mesh space dimension (2d or 3d space)
#if   __DIM == __SFIELD
#if   __MESH_DIM  == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_2d_sca_s(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_2d_sca_d(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_2d_sca_sc(fv,to_topo,from_mesh,to_mesh,  &
     &                                   ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_2d_sca_dc(fv,to_topo,from_mesh,to_mesh,  &
     &                                   ghostsize,maptype,info,mask)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_map_field_2d_sca_i(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_2d_sca_l(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#endif 
#elif __MESH_DIM  == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_3d_sca_s(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_3d_sca_d(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_3d_sca_sc(fv,to_topo,from_mesh,to_mesh,  &
     &                                   ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_3d_sca_dc(fv,to_topo,from_mesh,to_mesh,  &
     &                                   ghostsize,maptype,info,mask)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_map_field_3d_sca_i(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_3d_sca_l(fv,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#endif 
#endif

#elif __DIM == __VFIELD
#if   __MESH_DIM  == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_2d_vec_s(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_2d_vec_d(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_2d_vec_sc(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                   ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_2d_vec_dc(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                   ghostsize,maptype,info,mask)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_map_field_2d_vec_i(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_2d_vec_l(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#endif 
#elif __MESH_DIM  == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_3d_vec_s(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_3d_vec_d(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_3d_vec_sc(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_3d_vec_dc(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_map_field_3d_vec_i(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_3d_vec_l(fv,lda,to_topo,from_mesh,to_mesh,  &
     &                                  ghostsize,maptype,info,mask)
#endif 
#endif
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      USE ppm_module_map_field_send
      USE ppm_module_map_field_pop_2d
      USE ppm_module_map_field_pop_3d
      USE ppm_module_map_field_push_2d
      USE ppm_module_map_field_push_3d
      USE ppm_module_map_field_global
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __MESH_DIM  == __2D
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:)     , POINTER      :: fv
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:)     , POINTER      :: fv
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:)  , POINTER      :: fv
#else
      REAL(MK), DIMENSION(:,:,:)     , POINTER      :: fv
#endif
#elif __MESH_DIM  == __3D
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:)   , POINTER      :: fv
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:)   , POINTER      :: fv
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER      :: fv
#else
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER      :: fv
#endif
#endif

#elif __DIM == __VFIELD
#if   __MESH_DIM  == __2D
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:)     , POINTER    :: fv
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:)     , POINTER    :: fv
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:)  , POINTER    :: fv
#else
      REAL(MK), DIMENSION(:,:,:,:)     , POINTER    :: fv
#endif
#elif __MESH_DIM  == __3D
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:,:)   , POINTER    :: fv
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:,:)   , POINTER    :: fv
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:,:), POINTER    :: fv
#else
      REAL(MK), DIMENSION(:,:,:,:,:)   , POINTER    :: fv
#endif
#endif
#endif

#if   __MESH_DIM  == __2D
      LOGICAL, DIMENSION(:,:,:), POINTER, OPTIONAL  :: mask
#elif __MESH_DIM  == __3D
      LOGICAL, DIMENSION(:,:,:,:), POINTER,OPTIONAL :: mask
#endif

#if   __DIM == __VFIELD
      INTEGER                       , INTENT(IN   ) :: lda
#endif
      INTEGER                       , INTENT(IN   ) :: to_topo
      INTEGER                       , INTENT(IN   ) :: to_mesh
      INTEGER                       , INTENT(IN   ) :: from_mesh
      INTEGER                       , INTENT(IN   ) :: maptype
      INTEGER, DIMENSION(:)         , INTENT(IN   ) :: ghostsize
      INTEGER                       , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)              :: t0
      INTEGER               :: i
      CHARACTER(ppm_char)   :: mesg
      LOGICAL               :: valid
#if   __DIM == __SFIELD
      INTEGER, PARAMETER    :: lda = 1
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_field',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_field_topoid .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_no_topo,'ppm_map_field',  &
     &            'No field topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (to_topo .GT. 0) THEN
              CALL ppm_check_topoid(ppm_param_id_user,to_topo,valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &                'Invalid topology ID as to_topo specified!',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
#if   __DIM == __VFIELD
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (ghostsize(i) .LT. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &                'ghostsize must be >=0 in all dimensions',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine which mapping is required
      !-------------------------------------------------------------------------
      IF     (maptype.EQ.ppm_param_map_global) THEN
         !----------------------------------------------------------------------
         !  Full blast mapping; find the mesh points on other processors
         !----------------------------------------------------------------------
         ! if there is still some data left in the buffer, warn the user
         IF (ppm_buffer_set .GT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_map_incomp,'ppm_map_field',  &
     &           'Buffer was not empty. Possible loss of data!',__LINE__,info)
         ENDIF

         ! get new destination topology from user
         IF (to_topo .LE. 0) THEN
             ppm_target_topoid = ppm_field_topoid
         ELSE
             ppm_target_topoid = ppm_internal_topoid(to_topo)
         ENDIF

         !----------------------------------------------------------------------
         !  Check mesh id for validity
         !----------------------------------------------------------------------
         IF (ppm_debug .GT. 0) THEN
            IF (to_mesh .GT. 0) THEN
                CALL ppm_check_meshid(ppm_param_id_user,to_mesh,     &
     &              ppm_target_topoid,valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &                  'Destination mesh ID (to_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
            IF (from_mesh .GT. 0) THEN
                CALL ppm_check_meshid(ppm_param_id_user,from_mesh,     &
     &              ppm_field_topoid,valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &                  'Source mesh ID (from_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
         ENDIF

         ! get new source mesh from user
         IF (from_mesh .LE. 0) THEN
             ppm_source_meshid = 1
         ELSE
             ppm_source_meshid = ppm_meshid(ppm_field_topoid)%internal(from_mesh)
         ENDIF
         ! get new destination mesh from user
         IF (to_mesh .LE. 0) THEN
             ppm_target_meshid = 1
         ELSE
             ppm_target_meshid = ppm_meshid(ppm_target_topoid)%internal(to_mesh)
         ENDIF

         CALL ppm_map_field_global(ppm_target_topoid,ppm_source_meshid, &
    &                              ppm_target_meshid,info)
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_push) THEN
         !----------------------------------------------------------------------
         !  Add the data to the stack (push)
         !----------------------------------------------------------------------
         ! push can only be called for global
         IF (ppm_map_type .NE. ppm_param_map_global) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &           'push can only be called after _global',__LINE__,info)
             GOTO 9999
         ENDIF

         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0 .OR. ppm_target_meshid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_field',  &
     &           'Skipping push!',__LINE__,info)
             GOTO 9999
         ENDIF

         IF (PRESENT(mask)) THEN
#if   __MESH_DIM == __2D
             CALL ppm_map_field_push_2d(fv,lda,ppm_source_meshid,info,mask)
#elif __MESH_DIM == __3D
             CALL ppm_map_field_push_3d(fv,lda,ppm_source_meshid,info,mask)
#endif
         ELSE
#if   __MESH_DIM == __2D
             CALL ppm_map_field_push_2d(fv,lda,ppm_source_meshid,info)
#elif __MESH_DIM == __3D
             CALL ppm_map_field_push_3d(fv,lda,ppm_source_meshid,info)
#endif
         ENDIF
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_pop) THEN
         !----------------------------------------------------------------------
         !  Extract the data from the stack (pop) and replace current values
         !----------------------------------------------------------------------
         ! pop can only be called for global
         IF (ppm_map_type .NE. ppm_param_map_global) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &           'pop can only be called after _global',__LINE__,info)
             GOTO 9999
         ENDIF

         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0 .OR. ppm_target_meshid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_field',  &
     &           'Skipping pop!',__LINE__,info)
             GOTO 9999
         ENDIF

         ! skip if the buffer is empty
         IF (ppm_buffer_set .LT. 1) THEN
             IF (ppm_debug .GT. 1) THEN
                 CALL ppm_write(ppm_rank,'ppm_map_field',  &
     &               'Buffer is empty: skipping pop!',info)
             ENDIF
             GOTO 9999
         ENDIF

         IF (PRESENT(mask)) THEN
#if   __MESH_DIM == __2D
             CALL ppm_map_field_pop_2d(fv,lda,ppm_target_meshid,ghostsize,  &
     &           ppm_param_pop_replace,info,mask)
#elif __MESH_DIM == __3D
             CALL ppm_map_field_pop_3d(fv,lda,ppm_target_meshid,ghostsize,  &
     &           ppm_param_pop_replace,info,mask)
#endif
         ELSE
#if   __MESH_DIM == __2D
             CALL ppm_map_field_pop_2d(fv,lda,ppm_target_meshid,ghostsize,  &
     &           ppm_param_pop_replace,info)
#elif __MESH_DIM == __3D
             CALL ppm_map_field_pop_3d(fv,lda,ppm_target_meshid,ghostsize,  &
     &           ppm_param_pop_replace,info)
#endif
         ENDIF
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_send) THEN
         !----------------------------------------------------------------------
         !  Add the last data and send/recv the packages 
         !----------------------------------------------------------------------
         ! send can only be called for global
         IF (ppm_map_type .NE. ppm_param_map_global) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_map_field',  &
     &           'send can only be called after _global',__LINE__,info)
             GOTO 9999
         ENDIF

         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0 .OR. ppm_target_meshid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_field',  &
     &           'Skipping send!',__LINE__,info)
             GOTO 9999
         ENDIF

         ! skip if the buffer is empty
         IF (ppm_buffer_set .LT. 1) THEN
             IF (ppm_debug .GT. 1) THEN
                 CALL ppm_write(ppm_rank,'ppm_map_field',  &
     &               'Buffer is empty: skipping send!',info)
             ENDIF
             GOTO 9999
         ENDIF

         ! If mask was specified, push the mask as well
         IF (PRESENT(mask)) THEN
#if   __MESH_DIM == __2D
             CALL ppm_map_field_push_2d(mask,1,ppm_source_meshid,info)
#elif __MESH_DIM == __3D
             CALL ppm_map_field_push_3d(mask,1,ppm_source_meshid,info)
#endif
             IF (info .NE. 0) GOTO 9999
         ENDIF

         CALL ppm_map_field_send(info)
         IF (info.NE.0) GOTO 9999

         ! now the field is on the new topology
         ppm_field_topoid = ppm_target_topoid

         ! If mask was specified, pop the mask as well
         IF (PRESENT(mask)) THEN
#if   __MESH_DIM == __2D
             CALL ppm_map_field_pop_2d(mask,1,ppm_target_meshid,  &
     &           ghostsize,ppm_param_pop_replace,info)
#elif __MESH_DIM == __3D
             CALL ppm_map_field_pop_3d(mask,1,ppm_target_meshid,  &
     &           ghostsize,ppm_param_pop_replace,info)
#endif
             IF (info .NE. 0) GOTO 9999
         ENDIF
         
      ELSEIF (maptype.EQ.ppm_param_map_cancel) THEN
         !----------------------------------------------------------------------
         !  Cancel the mapping in progress and reset everything. No need to
         !  reset arrays since thay are alloc_fit-ed the next time a
         !  mapping is called. Their deallocation only happens in
         !  ppm_finalize anyway.
         !----------------------------------------------------------------------
         ppm_target_topoid = -1
         ppm_source_meshid = -1
         ppm_target_meshid = -1
         !----------------------------------------------------------------------
         !  These might actually not be needed since the mapping routines
         !  would in principle be expected to reset them. They are here for
         !  safety reasons.
         !----------------------------------------------------------------------
         ppm_buffer_set    = 0
         ppm_nsendbuffer   = 0
         ppm_nrecvbuffer   = 0

      ELSE
         !----------------------------------------------------------------------
         !  Unknow mapping
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_map_field',    &
     &       'Unknown mapping action specified',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field',t0,info)
      RETURN
#if   __DIM == __SFIELD
#if   __MESH_DIM  == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_2d_sca_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_2d_sca_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_2d_sca_dc
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_2d_sca_l
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_2d_sca_i
#endif
#elif __MESH_DIM  == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_3d_sca_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_3d_sca_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_3d_sca_dc
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_3d_sca_l
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_3d_sca_i
#endif
#endif

#elif __DIM == __VFIELD
#if   __MESH_DIM  == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_2d_vec_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_2d_vec_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_2d_vec_dc
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_2d_vec_l
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_2d_vec_i
#endif
#elif __MESH_DIM  == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_3d_vec_d
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_3d_vec_sc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_3d_vec_dc
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_3d_vec_l
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_3d_vec_i
#endif
#endif
#endif
