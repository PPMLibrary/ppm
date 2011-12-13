MODULE ppm_module_particles_typedef

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2


USE ppm_module_typedef

IMPLICIT NONE

TYPE pnt_array_1d_i
    INTEGER, DIMENSION(:), POINTER               :: vec => NULL()
    !!! array where the integer-valued property is stored
    CHARACTER(len=ppm_char)                      :: name
    !!! name of the integer-valued property
    LOGICAL                                      :: has_ghosts
    !!! true if ghost values are up-to-date
    LOGICAL                                      :: is_mapped
    !!! true if there is a one-to-one mapping with the particles
    LOGICAL                                      :: map_parts
    !!! true if partial mappings are desired for this property (default)
    LOGICAL                                      :: map_ghosts
    !!! true if ghost mappings are desired for this property (default)
END TYPE pnt_array_1d_i


#define  __KIND __SINGLE_PRECISION
#define  DTYPE(a) a/**/_s
#define  prec ppm_kind_single
#define  _prec _ppm_kind_single
#include "part/particles_typedef.inc"

#define  __KIND __DOUBLE_PRECISION
#define  DTYPE(a) a/**/_d
#define  prec ppm_kind_double
#define  _prec _ppm_kind_double
#include "part/particles_typedef.inc"

END MODULE ppm_module_particles_typedef

