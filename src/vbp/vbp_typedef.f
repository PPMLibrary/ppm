TYPE,EXTENDS(DTYPE(ppm_t_particles)) :: DTYPE(ppm_t_vbp)
    !!! an extension of the Particle set data structure
    !!! for Variable-Blob Particles

    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: rcp => NULL()
    !!! pointer to the property that stores the cutoff radius
        CONTAINS
    !        PRIVATE
           PROCEDURE :: create => DTYPE(vbp_create)
           PROCEDURE :: updated_cutoff   => DTYPE(vbp_updated_cutoff)
           PROCEDURE :: create_neighlist => DTYPE(vbp_neigh_create)
           PROCEDURE :: comp_neighlist   => DTYPE(vbp_neighlist)


END TYPE DTYPE(ppm_t_vbp)

#undef   MK
#undef   _MK
#undef   DTYPE
#undef   CTYPE
