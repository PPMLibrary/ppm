!!----------------------------------------------------------------------
!! Differential operators
!!----------------------------------------------------------------------
TYPE DTYPE(ppm_t_opdesc)
    !!! Data structure to describe differential operator
    !!!  It only contains semantic information on the operator
    INTEGER, DIMENSION(:), POINTER                 :: degree =>NULL()
    !!! degree of each term in the linear combination of differential ops 
    INTEGER, DIMENSION(:), POINTER                 :: order =>NULL()
    !!! approximation order of each term 
    REAL(MK), DIMENSION(:), POINTER                :: coeffs =>NULL()
    !!! array where the coefficients in linear combinations of 
    !!! differential ops are stored
    INTEGER                                        :: nterms
    !!! number of terms
    CHARACTER(LEN=ppm_char)                        :: name
    !!! name of the vector-valued property

    CONTAINS
    PROCEDURE :: create => DTYPE(desc_create)
    PROCEDURE :: destroy => DTYPE(desc_destroy)

END TYPE DTYPE(ppm_t_opdesc)

TYPE DTYPE(ppm_t_operator)
    !!! Data structure containing all diff operators for a particle set
    !!! 
    REAL(MK),DIMENSION(:,:),               POINTER :: ker => NULL()
    !!! where the operators are stored
    TYPE(DTYPE(ppm_t_opdesc)),             POINTER :: desc => NULL()
    !!! small matrices that describe what each operator does
    INTEGER                                        :: P_id = 0
    !!! Id of the set of particles that this operator takes data from.
    !!! The default, 0, stands for "self" (the operator is computed
    !!! on the same set of particles than the one which contains the data).
    INTEGER                                        :: neigh_id = 1
    !!! Id of the neighbour list that should be used
    !!! The default, 1, refers to "self": the list of neighbours within
    !!! the same set of particles. 

    LOGICAL, DIMENSION(ppm_param_length_opsflags)  :: flags
    !!! logical flags
    !!!    ppm_ops_inc_ghosts
    !!!           true if the operator should be computed for ghost 
    !!!           particles too.  Note that the resulting values 
    !!!           will be wrong for the ghost particles
    !!!           that have some neighbours outside the ghost layers. 
    !!!           Default is false.
    !!!    ppm_ops_interp
    !!!          true if the op interpolates data from one set of particles
    !!!    ppm_ops_iscomputed
    !!!          true if the operator has been computed and is uptodate
    !!!    ppm_ops_isdefined
    !!!          true if the operator has been defined 
    !!!    ppm_ops_vector
    !!!          true if each term represents a component (ie the result
    !!!          of the operator should be a vector field, like for gradients)
    !!!          false if the components are added up (like for the divergence)

    CONTAINS
    PROCEDURE :: create => DTYPE(op_create)
    PROCEDURE :: destroy => DTYPE(op_destroy)

END TYPE DTYPE(ppm_t_operator)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold operators
!----------------------------------------------------------------------
TYPE DTYPE(ppm_t_ptr_ops)
    TYPE(DTYPE(ppm_t_operator)), POINTER :: t => NULL()
END TYPE DTYPE(ppm_t_ptr_ops)

!----------------------------------------------------------------------
! Container for operators
!----------------------------------------------------------------------
#define CONTAINER DTYPE(ppm_c_operators)
#define __CONTAINER(a) DTYPE(ppm_c_operators)_/**/a
#define VEC_TYPE DTYPE(ppm_t_ptr_ops)
#define ITERATOR_TYPE DTYPE(ppm_t_operator)
#include "cont/extended_container_template.inc"

#undef   MK
#undef   _MK
#undef   DTYPE
