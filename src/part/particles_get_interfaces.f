#define WRAP(a) a

#define __FUNCNAME DTYPE(WRAP(DATANAME)_check)_
      SUBROUTINE __FUNCNAME(this,wp,info)
          IMPORT DTYPE(ppm_t_part_prop)_, ppm_kind_single
          IMPORT ppm_kind_double, ppm_kind_int64
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_part_prop)_)                                :: this
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(IN   ) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(IN   ) :: wp
#endif
          INTEGER,                                        INTENT(  OUT) :: info
      END SUBROUTINE
#undef __FUNCNAME
#undef __MYTYPE

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get_prop)_
      SUBROUTINE __FUNCNAME(this,discr_data,wp,info,with_ghosts,read_only,skip_checks)
          IMPORT DTYPE(ppm_t_particles)_, ppm_kind_single,ppm_kind_double
          IMPORT DTYPE(ppm_t_part_prop)_, ppm_kind_int64,ppm_t_field_
          IMPORT ppm_t_discr_data
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles)_)                                :: this
          !CLASS(DTYPE(ppm_t_part_prop)_),POINTER  :: discr_data
          CLASS(ppm_t_discr_data),                        INTENT(INOUT) :: discr_data
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          INTEGER,                                        INTENT(  OUT) :: info
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: with_ghosts
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: skip_checks
      END SUBROUTINE
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set_prop)_
      SUBROUTINE __FUNCNAME(this,discr_data,wp,info,read_only,ghosts_ok)
          IMPORT DTYPE(ppm_t_particles)_,ppm_kind_single,ppm_kind_double
          IMPORT ppm_t_discr_data, ppm_kind_int64,ppm_t_field_
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles)_)                                :: this
          !CLASS(DTYPE(ppm_t_part_prop)_),POINTER  :: discr_data
          CLASS(ppm_t_discr_data),                        INTENT(INOUT) :: discr_data
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          INTEGER,                                        INTENT(  OUT) :: info
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: ghosts_ok
      END SUBROUTINE
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get_field)_
      SUBROUTINE __FUNCNAME(this,Field,wp,info,with_ghosts,read_only,skip_checks)
          IMPORT DTYPE(ppm_t_particles)_, ppm_kind_single
          IMPORT ppm_kind_double,ppm_kind_int64,ppm_t_field_
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles)_)                                :: this
          CLASS(ppm_t_field_),                            INTENT(IN   ) :: Field
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          INTEGER,                                        INTENT(  OUT) :: info
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: with_ghosts
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: skip_checks
      END SUBROUTINE
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set_field)_
      SUBROUTINE __FUNCNAME(this,Field,wp,info,read_only,ghosts_ok)
          IMPORT DTYPE(ppm_t_particles)_,ppm_kind_single
          IMPORT ppm_kind_double,ppm_kind_int64,ppm_t_field_
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles)_)                                :: this
          CLASS(ppm_t_field_),                            INTENT(IN   ) :: Field
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          INTEGER,                                        INTENT(  OUT) :: info
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: ghosts_ok
      END SUBROUTINE
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get)_
#define __CHECKTYPE DTYPE(WRAP(DATANAME)_check)
      SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,info,with_ghosts,read_only,skip_checks)
          IMPORT DTYPE(ppm_t_particles)_, ppm_kind_single
          IMPORT ppm_kind_double, ppm_kind_int64
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles)_)                                :: Pc
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          INTEGER,                                        INTENT(IN   ) :: ppt_id
          INTEGER,                                        INTENT(  OUT) :: info
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: with_ghosts
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: skip_checks
      END SUBROUTINE
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set)_
      SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,info,read_only,ghosts_ok)
          IMPORT DTYPE(ppm_t_particles)_,ppm_kind_single
          IMPORT  ppm_kind_double, ppm_kind_int64
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles)_)                                :: Pc
#if   __DIM == 1
          __TYPE, DIMENSION(:),   POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          INTEGER,                                        INTENT(IN   ) :: ppt_id
          INTEGER,                                        INTENT(  OUT) :: info
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: ghosts_ok
      END SUBROUTINE
#undef __FUNCNAME

#undef DATANAME
#undef __TYPE
