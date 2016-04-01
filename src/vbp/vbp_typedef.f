      TYPE,EXTENDS(DTYPE(ppm_t_particles)) :: DTYPE(ppm_t_vbp)
          !!! an extension of the Particle set data structure
          !!! for Variable-Blob Particles

          LOGICAL                                 :: adaptive
          !!! true if the particles have their own cutoff radii
          !!! in this case, the cutoff will be stored in a property called rcp
          !!! Otherwise, the cutoff is a field of each neighbour list data object.
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: rcp => NULL()
          !!! pointer to the property that stores the cutoff radius

      CONTAINS

          PROCEDURE :: create             => DTYPE(vbp_create)
          PROCEDURE :: destroy            => DTYPE(vbp_destroy)
          PROCEDURE :: updated_cutoff     => DTYPE(vbp_updated_cutoff)
          PROCEDURE :: set_cutoff         => DTYPE(vbp_set_cutoff)
          PROCEDURE :: set_varying_cutoff => DTYPE(vbp_set_varying_cutoff)
          PROCEDURE :: create_neighlist   => DTYPE(vbp_neigh_create)
          PROCEDURE :: comp_neighlist     => DTYPE(vbp_neighlist)

      END TYPE DTYPE(ppm_t_vbp)

#undef   MK
#undef   _MK
#undef   DTYPE
#undef   CTYPE
