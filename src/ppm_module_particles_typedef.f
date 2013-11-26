      MODULE ppm_module_particles_typedef

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __REAL 3
#define __COMPLEX 4
#define __INTEGER 5
#define __LONGINT 6
#define __LOGICAL 7
#define __CHAR 8

      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_interfaces
      USE ppm_module_topo_typedef
      USE ppm_module_mapping_typedef

      IMPLICIT NONE

      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------

      INTEGER                        :: ppm_particles_seedsize
      INTEGER, DIMENSION(:), POINTER :: ppm_particles_seed => NULL()

#define  DTYPE(a) a/**/_s
#define  CTYPE(a) a/**/_sc
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "part/particles_typedef.f"

#define  DTYPE(a) a/**/_d
#define  CTYPE(a) a/**/_dc
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "part/particles_typedef.f"

      !----------------------------------------------------------------------
      ! DATA TYPE FOR COLLECTIVE STORAGE for particles
      !----------------------------------------------------------------------
      TYPE,EXTENDS(ppm_c_particles_s) :: ppm_vc_particles_s
      CONTAINS
          PROCEDURE :: vpush   => ppm_vc_particles_s_push
          PROCEDURE :: vremove => ppm_vc_particles_s_remove
      ENDTYPE
      TYPE,EXTENDS(ppm_c_particles_d) :: ppm_vc_particles_d
      CONTAINS
          PROCEDURE :: vpush   => ppm_vc_particles_d_push
          PROCEDURE :: vremove => ppm_vc_particles_d_remove
      ENDTYPE

      !----------------------------------------------------------------------
      ! DATA STORAGE for the particles
      !----------------------------------------------------------------------
      TYPE(ppm_vc_particles_s) :: ppm_part_s
      TYPE(ppm_vc_particles_d) :: ppm_part_d

      CHARACTER(LEN=ppm_char)        :: cbuf
      INTEGER, PRIVATE, DIMENSION(3) :: ldc
      !!! Number of elements in all dimensions for allocation

      CONTAINS

#include "part/ppm_particles_helpers.f"

#define DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_1d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_1d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_single)
#define __MYTYPE __REAL
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_single)
#define __MYTYPE __COMPLEX
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/particles_typeproc.f"
#include "part/ppm_part_neighlists_get.f"
#include "part/part_interp_to_mesh.f"
#undef  DEFINE_MK
# define DEFINE_MK() parameter(MK,<#INTEGER#>,ppm_kind_single)
#include "part/part_remesh.f"
#include "part/part_interp_to_mesh_all.f"
#include "part/particles_from_mesh.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND


#define DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define __DIM 1
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_1d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_1d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_double)
#define __MYTYPE __REAL
#define DATANAME data_1d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX
#define DATANAME data_1d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_1d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define __DIM 2
#define __TYPE INTEGER
#define __MYTYPE __INTEGER
#define DATANAME data_2d_i
#include "part/ppm_particles_get.f"
#define __TYPE INTEGER(ppm_kind_int64)
#define __MYTYPE __LONGINT
#define DATANAME data_2d_li
#include "part/ppm_particles_get.f"
#define __TYPE REAL(ppm_kind_double)
#define __MYTYPE __REAL
#define DATANAME data_2d_r
#include "part/ppm_particles_get.f"
#define __TYPE COMPLEX(ppm_kind_double)
#define __MYTYPE __COMPLEX
#define DATANAME data_2d_c
#include "part/ppm_particles_get.f"
#define __TYPE LOGICAL
#define __MYTYPE __LOGICAL
#define DATANAME data_2d_l
#include "part/ppm_particles_get.f"
#undef  __DIM

#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/particles_typeproc.f"
#include "part/ppm_part_neighlists_get.f"
#include "part/part_interp_to_mesh.f"
#undef  DEFINE_MK
# define DEFINE_MK() parameter(MK,<#INTEGER#>,ppm_kind_double)
#include "part/part_remesh.f"
#include "part/part_interp_to_mesh_all.f"
#include "part/particles_from_mesh.f"
#undef  DEFINE_MK

#undef  DTYPE
#undef  __KIND

      !PUSH
      SUBROUTINE ppm_vc_particles_s_push(this,element,info,id)
          !!! add an element into the collection

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_vc_particles_s), INTENT(INOUT) :: this
          CLASS(ppm_t_particles_s_), TARGET        :: element
          INTEGER,                   INTENT(  OUT) :: info
          INTEGER,         OPTIONAL, INTENT(  OUT) :: id
          !!! index of the element in the collection

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("ppm_vc_particles_s_push")

          !add the element at the end of the array
          this%min_id = 1
          this%max_id = this%max_id + 1
          this%nb = this%nb + 1
          IF (PRESENT(id)) id = this%max_id

          IF (this%max_id.GT.this%vec_size) THEN
             CALL this%grow_size(info)
             or_fail("could not grow ppm_vc_particles_s to a larger size")
          ENDIF

          IF (ASSOCIATED(this%vec(this%max_id)%t)) THEN
             fail("Pointer at position of new element is already associated. Something wrong in the Collection data structure")
          ENDIF

          this%vec(this%max_id)%t => element

          check_associated_noscope(<#this%vec(this%max_id)%t#>,"Pushing element into collection failed unexpectedly")

          end_subroutine()
      END SUBROUTINE ppm_vc_particles_s_push
      !REMOVE
      SUBROUTINE ppm_vc_particles_s_remove(this,info,element)
          !!! If element is present, remove it from the collection
          !!! else, remove the current element (as defined by the iterator pointer)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_vc_particles_s),         INTENT(INOUT) :: this
          INTEGER,                           INTENT(  OUT) :: info
          CLASS(ppm_t_particles_s_),OPTIONAL,INTENT(INOUT) :: element

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          INTEGER :: del_id
          INTEGER :: iter_id_save

          start_subroutine("ppm_vc_particles_s_remove")

          iter_id_save = this%iter_id

          IF (PRESENT(element)) THEN
             del_id = this%get_id(element)
             IF (del_id.LT.0) RETURN
          ELSE
             del_id = this%iter_id
          ENDIF

          !swap with the last non-empty element of the collection
          IF (this%max_id.GT.this%min_id) THEN
             this%vec(del_id)%t => this%vec(this%max_id)%t
             this%vec(this%max_id)%t => NULL()
          ELSE
             this%vec(del_id)%t => NULL()
          ENDIF

          this%nb = this%nb - 1
          this%max_id = this%max_id - 1
          this%iter_id = iter_id_save - 1
          IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0

          end_subroutine()
      END SUBROUTINE ppm_vc_particles_s_remove
      !PUSH
      SUBROUTINE ppm_vc_particles_d_push(this,element,info,id)
          !!! add an element into the collection

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_vc_particles_d), INTENT(INOUT) :: this
          CLASS(ppm_t_particles_d_), TARGET        :: element
          INTEGER,                   INTENT(  OUT) :: info
          INTEGER,         OPTIONAL, INTENT(  OUT) :: id
          !!! index of the element in the collection

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("ppm_vc_particles_d_push")

          !add the element at the end of the array
          this%min_id = 1
          this%max_id = this%max_id + 1
          this%nb = this%nb + 1
          IF (PRESENT(id)) id = this%max_id

          IF (this%max_id.GT.this%vec_size) THEN
             CALL this%grow_size(info)
             or_fail("could not grow ppm_vc_particles_d to a larger size")
          ENDIF

          IF (ASSOCIATED(this%vec(this%max_id)%t)) THEN
             fail("Pointer at position of new element is already associated. Something wrong in the Collection data structure")
          ENDIF

          this%vec(this%max_id)%t => element

          check_associated_noscope(<#this%vec(this%max_id)%t#>,"Pushing element into collection failed unexpectedly")

          end_subroutine()
      END SUBROUTINE ppm_vc_particles_d_push
      !REMOVE
      SUBROUTINE ppm_vc_particles_d_remove(this,info,element)
          !!! If element is present, remove it from the collection
          !!! else, remove the current element (as defined by the iterator pointer)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_vc_particles_d),         INTENT(INOUT) :: this
          INTEGER,                           INTENT(  OUT) :: info
          CLASS(ppm_t_particles_d_),OPTIONAL,INTENT(INOUT) :: element

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          INTEGER :: del_id
          INTEGER :: iter_id_save

          start_subroutine("ppm_vc_particles_d_remove")

          iter_id_save = this%iter_id

          IF (PRESENT(element)) THEN
             del_id = this%get_id(element)
             IF (del_id.LT.0) RETURN
          ELSE
             del_id = this%iter_id
          ENDIF

          !swap with the last non-empty element of the collection
          IF (this%max_id.GT.this%min_id) THEN
             this%vec(del_id)%t => this%vec(this%max_id)%t
             this%vec(this%max_id)%t => NULL()
          ELSE
             this%vec(del_id)%t => NULL()
          ENDIF

          this%nb = this%nb - 1
          this%max_id = this%max_id - 1
          this%iter_id = iter_id_save - 1
          IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0

          end_subroutine()
      END SUBROUTINE ppm_vc_particles_d_remove

      END MODULE ppm_module_particles_typedef

