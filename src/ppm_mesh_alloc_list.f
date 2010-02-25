      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_mesh_alloc_list
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Does the (re)allocation of arrays of type
      !                 ppm_type_mesh_list. It offers the same allocation
      !                 types as the ppm_alloc routines for regular arrays.
      !
      !  Input        : lda(:)    (I) new length of mesh list array.
      !                 iopt      (I) alloc action. One of:
      !                                 ppm_param_alloc_fit_preserve
      !                                 ppm_param_alloc_fit
      !                                 ppm_param_alloc_grow_preserve
      !                                 ppm_param_alloc_grow
      !                                 ppm_param_dealloc
      !
      !  Input/output : mesh_list (T) array of TYPE(ppm_type_mesh_list)
      !                               which is to be (re)allocated.
      !
      !  Output       : info      (I) Return status. 0 if everything OK.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mesh_alloc_list.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2004/10/01 16:33:37  ivos
      !  cosmetics.
      !
      !  Revision 1.6  2004/10/01 16:09:09  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/07/26 11:46:55  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.4  2004/07/26 07:42:47  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.3  2004/02/24 14:13:19  ivos
      !  Extended functionality: added ppm_param_alloc_fit, ppm_param_alloc_grow,
      !  ppm_param_alloc_grow_preserve and ppm_param_dealloc actions.
      !
      !  Revision 1.2  2004/02/19 18:01:47  ivos
      !  Now uses pointers to avoid explicit copying of data from arrays.
      !
      !  Revision 1.1  2004/02/19 14:20:26  ivos
      !  Added routine ppm_mesh_alloc for (re)allocation of mesh number list
      !  and mesh definition user-type arrays. The corresponding code has been
      !  removed from ppm_topo_mkfield and ppm_mesh_store and the Makefile
      !  updated. The new routines are in the ppm_module_mesh.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_alloc_list(mesh_list,lda,iopt,info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh, ONLY: ppm_type_mesh_list
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , DIMENSION(:), INTENT(IN   ) :: lda
      INTEGER                               , INTENT(IN   ) :: iopt
      TYPE(ppm_type_mesh_list), DIMENSION(:), POINTER :: mesh_list
      INTEGER                               , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER            :: i,ldc
      REAL(ppm_kind_double) :: t0
      TYPE(ppm_type_mesh_list), DIMENSION(:), POINTER :: work_list
      LOGICAL            :: lcopy,lalloc,lrealloc,ldealloc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_alloc_list',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0 .AND. iopt .NE. ppm_param_dealloc) THEN
          IF (SIZE(lda,1) .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_alloc_list',  &
     &            'lda must be at least of length 1',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_alloc_list',  &
     &            'lda(1) must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check allocation type
      !-------------------------------------------------------------------------
      lcopy    = .FALSE.
      lalloc   = .FALSE.
      lrealloc = .FALSE.
      ldealloc = .FALSE.
      IF (iopt .EQ. ppm_param_alloc_fit_preserve) THEN
          !---------------------------------------------------------------------
          !  Fit memory and preserve the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(mesh_list)) THEN
              ldc = SIZE(mesh_list)
              IF (ldc .NE. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  lcopy    = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_alloc_fit) THEN
          !---------------------------------------------------------------------
          !  Fit memory and discard the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(mesh_list)) THEN
              ldc = SIZE(mesh_list)
              IF (ldc .NE. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  ldealloc = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_alloc_grow_preserve) THEN
          !---------------------------------------------------------------------
          !  Fit memory and preserve the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(mesh_list)) THEN
              ldc = SIZE(mesh_list)
              IF (ldc .LT. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  lcopy    = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_alloc_grow) THEN
          !---------------------------------------------------------------------
          !  Fit memory and discard the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(mesh_list)) THEN
              ldc = SIZE(mesh_list)
              IF (ldc .LT. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  ldealloc = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_dealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate
          !---------------------------------------------------------------------
          IF (ASSOCIATED(mesh_list)) THEN
              ldc = SIZE(mesh_list)
              lrealloc = .TRUE.
              ldealloc = .TRUE.
          ENDIF
      ENDIF
              
      !-------------------------------------------------------------------------
      !  Perform the actual alloc action
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
          !---------------------------------------------------------------------
          !  Allocate new array with new size and nullify all members
          !---------------------------------------------------------------------
          ALLOCATE(work_list(lda(1)),STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_mesh_alloc_list',  &
     &            'new mesh list WORK_LIST',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,lda(1)
              NULLIFY(work_list(i)%internal)
              NULLIFY(work_list(i)%user)
          ENDDO
      ENDIF

      IF (lcopy) THEN
          !---------------------------------------------------------------------
          !  Save the old contents
          !---------------------------------------------------------------------
          DO i=1,MIN(ldc,lda(1))
              work_list(i)%user => mesh_list(i)%user
              work_list(i)%internal => mesh_list(i)%internal
          ENDDO
      ENDIF

      IF (ldealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate the old contents
          !---------------------------------------------------------------------
          DO i=1,ldc
              IF (ASSOCIATED(mesh_list(i)%internal)) THEN
                  DEALLOCATE(mesh_list(i)%internal,STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_list',  &
     &                  'internal number list MESH_LIST%INTERNAL',__LINE__,info)
                        GOTO 9999
                  ENDIF
                  NULLIFY(mesh_list(i)%internal)
              ENDIF
              IF (ASSOCIATED(mesh_list(i)%user)) THEN
                  DEALLOCATE(mesh_list(i)%user,STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_list',  &
     &                   'user number list MESH_LIST%USER',__LINE__,info)
                      GOTO 9999
                  ENDIF
                  NULLIFY(mesh_list(i)%user)
              ENDIF
          ENDDO
      ENDIF

      IF (lrealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate old pointer array
          !---------------------------------------------------------------------
          DEALLOCATE(mesh_list,STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_list',  &
     &            'mesh number list MESH_LIST',__LINE__,info)
              GOTO 9999
          ENDIF
          NULLIFY(mesh_list)
      ENDIF

      IF (lalloc) THEN
          !---------------------------------------------------------------------
          !  Point result to new array
          !---------------------------------------------------------------------
          mesh_list => work_list
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_alloc_list',t0,info)
      RETURN
      END SUBROUTINE ppm_mesh_alloc_list
