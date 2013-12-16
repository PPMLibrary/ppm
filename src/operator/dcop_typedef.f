      TYPE,EXTENDS(ppm_t_operator_discr) :: DTYPE(ppm_t_dcop)
          !!! Data structure containing all diff operators for a particle set
          !!!
          REAL(MK),DIMENSION(:,:),        POINTER :: ker => NULL()
          !!! where the operators are stored
          !!! small matrices that describe what each operator does
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: nn_sq => NULL()
          !!! Pointer to the property containing the squared nearest-neighbour
          !!! distance between particles. Needed only for interpolating operators.

      CONTAINS

          PROCEDURE :: create          => DTYPE(dcop_create)
          PROCEDURE :: destroy         => DTYPE(dcop_destroy)
          PROCEDURE :: compute         => DTYPE(dcop_compute)
          PROCEDURE :: comp_weights_2d => DTYPE(dcop_comp_weights_2d)
          PROCEDURE :: comp_weights_3d => DTYPE(dcop_comp_weights_3d)
          PROCEDURE :: comp_weights    => DTYPE(dcop_comp_weights)

      END TYPE

#undef DTYPE
#undef MK
#undef _MK
