      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_mesh_alloc_equi
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_alloc_equi(equi_mesh,lda,iopt,info)
      !!! Does the (re)allocation of arrays of type `ppm_type_equi_mesh`.
      !!!
      !!! It offers the same allocation types as ppm_alloc for regular arrays.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_typedef, ONLY: ppm_t_equi_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , DIMENSION(:  ), INTENT(IN   ) :: lda
      !!! New size of mesh definition array
      INTEGER                                 , INTENT(IN   ) :: iopt
      !!! Alloc action. One of:
      !!!
      !!! * ppm_param_alloc_fit_preserve
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_alloc_grow_preserve
      !!! * ppm_param_alloc_grow
      !!! * ppm_param_dealloc
      TYPE(ppm_t_equi_mesh)   , DIMENSION(:  ), POINTER :: equi_mesh
      !!! Array of TYPE(ppm_t_equi_mesh) which is to be (re)allocated.
      INTEGER                                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER            :: i,j
      INTEGER, DIMENSION(2) :: ldc
      REAL(ppm_kind_double) :: t0
      TYPE(ppm_t_equi_mesh)   , DIMENSION(:  ), POINTER :: work_mesh
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
        CALL check
        IF (info .NE. 0) GOTO 9999
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
              IF (ldc(1) .NE. lda(1)) THEN
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
              IF (ldc(1) .NE. lda(1)) THEN
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
              IF (ldc(1) .LT. lda(1)) THEN
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
              IF (ldc(1) .LT. lda(1)) THEN
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
          ALLOCATE(work_mesh(lda(1)),STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_mesh_alloc_equi',  &
     &            'new mesh WORK_MESH',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,lda(1)
              NULLIFY(work_mesh(i)%istart)
              NULLIFY(work_mesh(i)%nnodes)
              NULLIFY(work_mesh(i)%Nm)
              work_mesh(i)%ID = 0
          ENDDO
      ENDIF

      IF (lcopy) THEN
          !---------------------------------------------------------------------
          !  Save the old contents
          !---------------------------------------------------------------------
          DO i=1,MIN(ldc(1),lda(1))
              work_mesh(i)%istart => equi_mesh(i)%istart
              work_mesh(i)%nnodes => equi_mesh(i)%nnodes
              work_mesh(i)%Nm => equi_mesh(i)%Nm
              work_mesh(i)%ID = equi_mesh(i)%ID
          ENDDO
      ENDIF

      IF (ldealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate the old contents
          !---------------------------------------------------------------------
          DO i=1,ldc(1)
              IF (ASSOCIATED(equi_mesh(i)%istart)) THEN
                  DEALLOCATE(equi_mesh(i)%istart,STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',&
     &                   'mesh start EQUI_MESH%ISTART',__LINE__,info)
                  ENDIF
                  NULLIFY(equi_mesh(i)%istart)
              ENDIF
              IF (ASSOCIATED(equi_mesh(i)%nnodes)) THEN
                  DEALLOCATE(equi_mesh(i)%nnodes,STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',&
     &                   'mesh size EQUI_MESH%NNODES',__LINE__,info)
                  ENDIF
                  NULLIFY(equi_mesh(i)%nnodes)
              ENDIF
              IF (ASSOCIATED(equi_mesh(i)%Nm)) THEN
                  DEALLOCATE(equi_mesh(i)%Nm,STAT=info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_dealloc,'ppm_mesh_alloc_equi',&
     &                   'global mesh size EQUI_MESH%NM',__LINE__,info)
                  ENDIF
                  NULLIFY(equi_mesh(i)%Nm)
              ENDIF
              equi_mesh(i)%ID = 0
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
      CONTAINS
      SUBROUTINE check
          IF (lda(1) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_alloc_equi',  &
     &            'lda(1) must be >= 0',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_alloc_equi
