      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_mesh_alloc_equi
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Does the (re)allocation of arrays of type
      !                 ppm_type_equi_mesh. It offers the same allocation
      !                 types as ppm_alloc for regular arrays.
      !
      !  Input        : lda(:)    (I) new size of mesh definition array in
      !                               both dimensions
      !                 iopt      (I) alloc action. One of:
      !                                 ppm_param_alloc_fit_preserve
      !                                 ppm_param_alloc_fit
      !                                 ppm_param_alloc_grow_preserve
      !                                 ppm_param_alloc_grow
      !                                 ppm_param_dealloc
      !
      !  Input/output : equi_mesh (T) array of TYPE(ppm_type_equi_mesh)
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
      !  $Log: ppm_mesh_alloc_equi.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2004/12/03 17:36:16  michaebe
      !  Fixed a bug where the topoid wasnt copied in fit_preserve
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
      !  Revision 1.1  2004/02/19 14:20:27  ivos
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

      SUBROUTINE ppm_mesh_alloc_equi(equi_mesh,lda,iopt,info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Module
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh, ONLY: ppm_type_equi_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , DIMENSION(:  ), INTENT(IN   ) :: lda
      INTEGER                                 , INTENT(IN   ) :: iopt
      TYPE(ppm_type_equi_mesh), DIMENSION(:,:), POINTER :: equi_mesh
      INTEGER                                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER            :: i,j
      INTEGER, DIMENSION(2) :: ldc
      REAL(ppm_kind_double) :: t0
      TYPE(ppm_type_equi_mesh), DIMENSION(:,:), POINTER :: work_mesh
      LOGICAL            :: lcopy,lalloc,lrealloc,ldealloc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_alloc_equi',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0 .AND. iopt .NE. ppm_param_dealloc) THEN
          IF (SIZE(lda,1) .LT. 2) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_alloc_equi',  &
     &            'lda must be at least of length 2',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_alloc_equi',  &
     &            'lda(1) must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_alloc_equi',  &
     &            'lda(2) must be >= 0',__LINE__,info)
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
          IF (ASSOCIATED(equi_mesh)) THEN
              ldc(1) = SIZE(equi_mesh,1)
              ldc(2) = SIZE(equi_mesh,2)
              IF ((ldc(1) .NE. lda(1)) .OR. (ldc(2) .NE. lda(2))) THEN
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
          IF (ASSOCIATED(equi_mesh)) THEN
              ldc(1) = SIZE(equi_mesh,1)
              ldc(2) = SIZE(equi_mesh,2)
              IF ((ldc(1) .NE. lda(1)) .OR. (ldc(2) .NE. lda(2))) THEN
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
          IF (ASSOCIATED(equi_mesh)) THEN
              ldc(1) = SIZE(equi_mesh,1)
              ldc(2) = SIZE(equi_mesh,2)
              IF ((ldc(1) .LT. lda(1)) .OR. (ldc(2) .LT. lda(2))) THEN
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
          IF (ASSOCIATED(equi_mesh)) THEN
              ldc(1) = SIZE(equi_mesh,1)
              ldc(2) = SIZE(equi_mesh,2)
              IF ((ldc(1) .LT. lda(1)) .OR. (ldc(2) .LT. lda(2))) THEN
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
          IF (ASSOCIATED(equi_mesh)) THEN
              ldc(1) = SIZE(equi_mesh,1)
              ldc(2) = SIZE(equi_mesh,2)
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
          ALLOCATE(work_mesh(lda(1),lda(2)),STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_mesh_alloc_equi',  &
     &            'new mesh WORK_MESH',__LINE__,info)
              GOTO 9999
          ENDIF
          DO j=1,lda(2)
              DO i=1,lda(1)
                  NULLIFY(work_mesh(i,j)%istart)
                  NULLIFY(work_mesh(i,j)%nnodes)
                  NULLIFY(work_mesh(i,j)%Nm)
              ENDDO
          ENDDO
      ENDIF

      IF (lcopy) THEN
          !---------------------------------------------------------------------
          !  Save the old contents
          !---------------------------------------------------------------------
          DO j=1,MIN(ldc(2),lda(2))
              DO i=1,MIN(ldc(1),lda(1))
                  work_mesh(i,j)%istart => equi_mesh(i,j)%istart
                  work_mesh(i,j)%nnodes => equi_mesh(i,j)%nnodes
                  work_mesh(i,j)%Nm => equi_mesh(i,j)%Nm
                  work_mesh(i,j)%topoid = equi_mesh(i,j)%topoid
              ENDDO
          ENDDO
      ENDIF

      IF (ldealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate the old contents
          !---------------------------------------------------------------------
          DO j=1,ldc(2)
              DO i=1,ldc(1)
                  IF (ASSOCIATED(equi_mesh(i,j)%istart)) THEN
                      DEALLOCATE(equi_mesh(i,j)%istart,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',&
     &                       'mesh start EQUI_MESH%ISTART',__LINE__,info)
                      ENDIF
                      NULLIFY(equi_mesh(i,j)%istart)
                  ENDIF
                  IF (ASSOCIATED(equi_mesh(i,j)%nnodes)) THEN
                      DEALLOCATE(equi_mesh(i,j)%nnodes,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',&
     &                       'mesh size EQUI_MESH%NNODES',__LINE__,info)
                      ENDIF
                      NULLIFY(equi_mesh(i,j)%nnodes)
                  ENDIF
                  IF (ASSOCIATED(equi_mesh(i,j)%Nm)) THEN
                      DEALLOCATE(equi_mesh(i,j)%Nm,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',&
     &                       'global mesh size EQUI_MESH%NM',__LINE__,info)
                      ENDIF
                      NULLIFY(equi_mesh(i,j)%Nm)
                  ENDIF
              ENDDO
          ENDDO
      ENDIF

      IF (lrealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate old pointer array
          !---------------------------------------------------------------------
          DEALLOCATE(equi_mesh,STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',  &
     &            'mesh data EQUI_MESH',__LINE__,info)
              GOTO 9999
          ENDIF
          NULLIFY(equi_mesh)
      ENDIF

      IF (lalloc) THEN
          !---------------------------------------------------------------------
          !  Point result to new array
          !---------------------------------------------------------------------
          equi_mesh => work_mesh
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_alloc_equi',t0,info)
      RETURN
      END SUBROUTINE ppm_mesh_alloc_equi
