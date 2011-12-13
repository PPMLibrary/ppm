MODULE ppm_module_particles

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __crash_on_null_pointers 1

USE ppm_module_particles_typedef
USE ppm_module_typedef
USE ppm_module_alloc
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_error
USE ppm_module_write
USE ppm_module_data, ONLY: ppm_dim,ppm_rank

IMPLICIT NONE

    !----------------------------------------------------------------------
    ! Global variables and parameters
    !----------------------------------------------------------------------
    INTEGER, PARAMETER :: ppm_param_part_init_cartesian = 1
    INTEGER, PARAMETER :: ppm_param_part_init_random = 2

    INTEGER                               :: ppm_particles_seedsize
    INTEGER,  DIMENSION(:  ), POINTER     :: ppm_particles_seed => NULL()

    !----------------------------------------------------------------------
    ! Private variables for the module
    !----------------------------------------------------------------------
    !for debugging
    LOGICAL,PRIVATE                   :: verbose = .FALSE.

    CHARACTER(LEN=32)      :: line_of_stars='********************************'
    INTEGER, PRIVATE, DIMENSION(3)    :: ldc
    !!! Number of elements in all dimensions for allocation
    CHARACTER(LEN=ppm_char),PRIVATE   :: cbuf
    !!! buffer for writeouts

    !REAL(prec),DIMENSION(:),POINTER    :: tmp_cutoff


#define interface_s_d(a) \
    INTERFACE a ;\
        MODULE PROCEDURE a/**/_s; \
        MODULE PROCEDURE a/**/_d; \
    END INTERFACE

    
    interface_s_d(get_xp)
    interface_s_d(set_xp)
    interface_s_d(get_wpi)
    interface_s_d(set_wpi)
    interface_s_d(get_wps)
    interface_s_d(set_wps)
    interface_s_d(get_wpv)
    interface_s_d(set_wpv)
    interface_s_d(get_dcop)
    interface_s_d(set_dcop)
    interface_s_d(ppm_alloc_particles)
    interface_s_d(particles_dcop_deallocate)
    interface_s_d(particles_allocate_wpi)
    interface_s_d(particles_allocate_wps)
    interface_s_d(particles_allocate_wpv)
    interface_s_d(particles_mapping_global)
    interface_s_d(particles_mapping_partial)
    interface_s_d(particles_mapping_ghosts)
    interface_s_d(particles_apply_bc)
    interface_s_d(particles_move)
    interface_s_d(particles_have_moved)
    interface_s_d(particles_neighlists)
    interface_s_d(particles_neighlists_xset)
    interface_s_d(particles_updated_positions)
    interface_s_d(particles_updated_nb_part)
    interface_s_d(particles_update_cutoff)
    interface_s_d(particles_updated_cutoff)
    interface_s_d(particles_compute_hmin)
    interface_s_d(particles_initialize)
    interface_s_d(particles_io_xyz)
    interface_s_d(particles_dcop_define)
    interface_s_d(particles_dcop_free)
    interface_s_d(particles_dcop_apply)
    interface_s_d(particles_apply_dcops)
    interface_s_d(particles_print_stats)
    interface_s_d(particles_check_arrays)
    interface_s_d(particles_add)
    interface_s_d(particles_initialize2d)
    interface_s_d(particles_initialize3d)

CONTAINS

#include "part/ppm_particles_helpers.f"

#define  __KIND __SINGLE_PRECISION
#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "part/ppm_particles.f"

#define  __KIND __DOUBLE_PRECISION
#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "part/ppm_particles.f"


END MODULE ppm_module_particles
