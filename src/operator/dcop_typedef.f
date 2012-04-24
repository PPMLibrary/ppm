TYPE,EXTENDS(ppm_t_operator_discr) :: DTYPE(ppm_t_dcop)
    !!! Data structure containing all diff operators for a particle set
    !!! 
    REAL(MK),DIMENSION(:,:),               POINTER :: ker => NULL()
    !!! where the operators are stored
    !!! small matrices that describe what each operator does
    CONTAINS
    PROCEDURE create  => DTYPE(dcop_create)
    PROCEDURE destroy => DTYPE(dcop_destroy)
    PROCEDURE compute => DTYPE(dcop_compute)

END TYPE

#undef DTYPE
#undef MK
#undef _MK
