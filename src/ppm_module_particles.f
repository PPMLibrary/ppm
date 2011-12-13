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


    INTERFACE get_xp
        MODULE PROCEDURE get_xp_s
        MODULE PROCEDURE get_xp_d
    END INTERFACE
    INTERFACE set_xp
        MODULE PROCEDURE set_xp_s
        MODULE PROCEDURE set_xp_d
    END INTERFACE
    INTERFACE get_wpi
        MODULE PROCEDURE get_wpi_s
        MODULE PROCEDURE get_wpi_d
    END INTERFACE
    INTERFACE set_wpi
        MODULE PROCEDURE set_wpi_s
        MODULE PROCEDURE set_wpi_d
    END INTERFACE
    INTERFACE get_wps
        MODULE PROCEDURE get_wps_s
        MODULE PROCEDURE get_wps_d
    END INTERFACE
    INTERFACE set_wps
        MODULE PROCEDURE set_wps_s
        MODULE PROCEDURE set_wps_d
    END INTERFACE
    INTERFACE get_wpv
        MODULE PROCEDURE get_wpv_s
        MODULE PROCEDURE get_wpv_d
    END INTERFACE
    INTERFACE set_wpv
        MODULE PROCEDURE set_wpv_s
        MODULE PROCEDURE set_wpv_d
    END INTERFACE
    INTERFACE get_dcop
        MODULE PROCEDURE get_dcop_s
        MODULE PROCEDURE get_dcop_d
    END INTERFACE
    INTERFACE set_dcop
        MODULE PROCEDURE set_dcop_s
        MODULE PROCEDURE set_dcop_d
    END INTERFACE
    INTERFACE ppm_alloc_particles
        MODULE PROCEDURE ppm_alloc_particles_s
        MODULE PROCEDURE ppm_alloc_particles_d
    END INTERFACE
    INTERFACE particles_dcop_deallocate
        MODULE PROCEDURE particles_dcop_deallocate_s
        MODULE PROCEDURE particles_dcop_deallocate_d
    END INTERFACE
    INTERFACE particles_allocate_wpi
        MODULE PROCEDURE particles_allocate_wpi_s
        MODULE PROCEDURE particles_allocate_wpi_d
    END INTERFACE
    INTERFACE particles_allocate_wps
        MODULE PROCEDURE particles_allocate_wps_s
        MODULE PROCEDURE particles_allocate_wps_d
    END INTERFACE
    INTERFACE particles_allocate_wpv
        MODULE PROCEDURE particles_allocate_wpv_s
        MODULE PROCEDURE particles_allocate_wpv_d
    END INTERFACE
    INTERFACE particles_mapping_global
        MODULE PROCEDURE particles_mapping_global_s
        MODULE PROCEDURE particles_mapping_global_d
    END INTERFACE
    INTERFACE particles_mapping_partial
        MODULE PROCEDURE particles_mapping_partial_s
        MODULE PROCEDURE particles_mapping_partial_d
    END INTERFACE
    INTERFACE particles_mapping_ghosts
        MODULE PROCEDURE particles_mapping_ghosts_s
        MODULE PROCEDURE particles_mapping_ghosts_d
    END INTERFACE
    INTERFACE particles_apply_bc
        MODULE PROCEDURE particles_apply_bc_s
        MODULE PROCEDURE particles_apply_bc_d
    END INTERFACE
    INTERFACE particles_move
        MODULE PROCEDURE particles_move_s
        MODULE PROCEDURE particles_move_d
    END INTERFACE
    INTERFACE particles_have_moved
        MODULE PROCEDURE particles_have_moved_s
        MODULE PROCEDURE particles_have_moved_d
    END INTERFACE
    INTERFACE particles_neighlists
        MODULE PROCEDURE particles_neighlists_s
        MODULE PROCEDURE particles_neighlists_d
    END INTERFACE
    INTERFACE particles_neighlists_xset
        MODULE PROCEDURE particles_neighlists_xset_s
        MODULE PROCEDURE particles_neighlists_xset_d
    END INTERFACE
    INTERFACE particles_updated_positions
        MODULE PROCEDURE particles_updated_positions_s
        MODULE PROCEDURE particles_updated_positions_d
    END INTERFACE
    INTERFACE particles_updated_nb_part
        MODULE PROCEDURE particles_updated_nb_part_s
        MODULE PROCEDURE particles_updated_nb_part_d
    END INTERFACE
    INTERFACE particles_update_cutoff
        MODULE PROCEDURE particles_update_cutoff_s
        MODULE PROCEDURE particles_update_cutoff_d
    END INTERFACE
    INTERFACE particles_updated_cutoff
        MODULE PROCEDURE particles_updated_cutoff_s
        MODULE PROCEDURE particles_updated_cutoff_d
    END INTERFACE
    INTERFACE particles_compute_hmin
        MODULE PROCEDURE particles_compute_hmin_s
        MODULE PROCEDURE particles_compute_hmin_d
    END INTERFACE
    INTERFACE particles_initialize
        MODULE PROCEDURE particles_initialize_s
        MODULE PROCEDURE particles_initialize_d
    END INTERFACE
    INTERFACE particles_io_xyz
        MODULE PROCEDURE particles_io_xyz_s
        MODULE PROCEDURE particles_io_xyz_d
    END INTERFACE
    INTERFACE particles_dcop_define
        MODULE PROCEDURE particles_dcop_define_s
        MODULE PROCEDURE particles_dcop_define_d
    END INTERFACE
    INTERFACE particles_dcop_free
        MODULE PROCEDURE particles_dcop_free_s
        MODULE PROCEDURE particles_dcop_free_d
    END INTERFACE
    INTERFACE particles_dcop_apply
        MODULE PROCEDURE particles_dcop_apply_s
        MODULE PROCEDURE particles_dcop_apply_d
    END INTERFACE
    INTERFACE particles_apply_dcops
        MODULE PROCEDURE particles_apply_dcops_s
        MODULE PROCEDURE particles_apply_dcops_d
    END INTERFACE
    INTERFACE particles_print_stats
        MODULE PROCEDURE particles_print_stats_s
        MODULE PROCEDURE particles_print_stats_d
    END INTERFACE
    INTERFACE particles_check_arrays
        MODULE PROCEDURE particles_check_arrays_s
        MODULE PROCEDURE particles_check_arrays_d
    END INTERFACE
    INTERFACE particles_add
        MODULE PROCEDURE particles_add_s
        MODULE PROCEDURE particles_add_d
    END INTERFACE
    INTERFACE particles_initialize2d
        MODULE PROCEDURE particles_initialize2d_s
        MODULE PROCEDURE particles_initialize2d_d
    END INTERFACE
    INTERFACE particles_initialize3d
        MODULE PROCEDURE particles_initialize3d_s
        MODULE PROCEDURE particles_initialize3d_d
    END INTERFACE

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
