SUBROUTINE ppm_part_dcops(Particles,eta_id,c,info,&
        order_deriv,order_approx,islaplacian,iscartesian,&
        isinterp,Particles_old,nn_sq_id)

    USE ppm_module_data, ONLY: ppm_dim,ppm_rank
    USE ppm_module_error
    USE ppm_module_particles
    IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION 
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    TYPE(ppm_t_particles),   POINTER,    INTENT(INOUT)   :: Particles
    !!! particles
    INTEGER,                             INTENT(INOUT)   :: eta_id
    !!! id of the operator kernel
    REAL(MK),                            INTENT(IN   )   :: c
    !!! ratio h/epsilon
    INTEGER,                             INTENT(  OUT)   :: info
    !!! non-zero on output if some error occurred
    !---------------------------------------------------------
    ! Optional arguments
    !---------------------------------------------------------
    INTEGER,                   OPTIONAL, INTENT(IN   )   :: order_approx
    !!! order of approximation
    INTEGER,DIMENSION(ppm_dim),OPTIONAL, INTENT(IN   )   :: order_deriv
    !!! degree of the derivative that we want to approximate
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: islaplacian
    !!! special treatment for the laplacian operator. Default is FALSE
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: iscartesian
    !!! special treatment if particles are on a grid. Default is FALSE
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: isinterp
    !!! special treatment if the operator is interpolating. Default is FALSE
    TYPE(ppm_t_particles),OPTIONAL,POINTER,INTENT(IN   ) :: Particles_old
    !!! Second set of particles. If present, then the kernel uses data
    !!! on the old set (used to compute interpolation kernels e.g.)
    INTEGER,                   OPTIONAL, INTENT(IN   )   :: nn_sq_id
    !!! id where the nearest-neighbour distances within Particles_old 
    !!! are stored (they must have been already computed)

    info = 0

    IF (ppm_dim .EQ. 2) THEN
        CALL ppm_part_dcops_2d(Particles,eta_id,c,info,&
            order_deriv,order_approx,islaplacian,iscartesian,&
            isinterp,Particles_old,nn_sq_id)
    ELSE
        CALL ppm_part_dcops_3d(Particles,eta_id,c,info,&
            order_deriv,order_approx,islaplacian,iscartesian,&
            isinterp,Particles_old,nn_sq_id)
    ENDIF

    9999 CONTINUE ! jump here upon error

END SUBROUTINE ppm_part_dcops

#undef __KIND
