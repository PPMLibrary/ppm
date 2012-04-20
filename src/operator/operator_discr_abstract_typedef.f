
!!----------------------------------------------------------------------
!! DC operators
!!----------------------------------------------------------------------

!TYPE,ABSTRACT,EXTENDS(ppm_t_operator_discr_) :: DTYPE(ppm_t_dcop)_
    !!!! Data structure containing all diff operators for a particle set
    !!!! 
    !REAL(MK),DIMENSION(:,:),               POINTER :: ker => NULL()
    !!!! where the operators are stored
    !CLASS(DTYPE(ppm_t_opdesc)_),           POINTER :: desc => NULL()
    !!!! small matrices that describe what each operator does
    !CONTAINS
    !PROCEDURE(DTYPE(dcop_create)_),  DEFERRED :: create
    !PROCEDURE(DTYPE(dcop_destroy)_), DEFERRED :: destroy
    !PROCEDURE(DTYPE(dcop_compute)_), DEFERRED :: compute

!END TYPE DTYPE(ppm_t_dcop)_
!minclude define_abstract_collection_type(DTYPE(ppm_t_dcop)_)
