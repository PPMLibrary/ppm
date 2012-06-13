TYPE,EXTENDS(DTYPE(ppm_t_vbp)) :: DTYPE(ppm_t_sop)
    !!! Extension of the Variable-Blob Particle set data structure
    !!! for Self-Organizing Particles

    ! Adaptive particles
    LOGICAL                                         :: adaptive
    !!! true if the particles have their own cutoff radii
    !!! in this case, the cutoff will be stored in wps(rcp_id)%vec
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: D => NULL()
    !!! pointer to the property that stores D
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: Dtilde => NULL()
    !!! pointer to the property that stores D_tilde
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: adapt_wp => NULL()
    !!! pointer to the property on which the adaptation is based
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: adapt_wp_grad => NULL()
    !!! pointer to the gradient of adapt_wp
    LOGICAL                                         :: level_set
    !!! true if particles carry a level-set function
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: level => NULL()
    !!! index of the wps array where the level-set is stored
    !    INTEGER                                         :: level_old_id
    !!! index of the wps array where the level-set is backed up before adapt

    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: level_grad => NULL()
    !!! index of the wps array where the gradient of the level-set is stored
    !    INTEGER                                         :: level_grad_old_id
    !!! index of the wps array where the gradient of the level-set 
    !!! is backed up before adapt


    ! List of IDs of other adaptive Particle sets
    !TYPE(idList)                                    :: set_aPc


    ! Anisotropic particles
    LOGICAL                                         :: anisotropic
    !!! true if the particles have their own cutoff radii
    !!! in this case, the G tensor will be stored in wpv(G_id)%vec
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: G => NULL()
    !!! index where G is stored

        CONTAINS
    !        PRIVATE
           PROCEDURE :: DTYPE(sop_adapt_particles)


END TYPE DTYPE(ppm_t_sop)

TYPE DTYPE(sop_t_opts)
    !!! derived type for optional arguments to the adaptivity routine

    LOGICAL           :: level_set
    REAL(MK)          :: param_nb
    REAL(MK)          :: nb_width
    REAL(MK)          :: nb_width2
    REAL(MK)          :: nb_width_kill
    INTEGER           :: order_approx
    !!! order of approximation for the interpolation kernels
    REAL(MK)          :: c
    LOGICAL           :: check_dcops

    LOGICAL           :: write_pdb    ! writeout pdb files
    LOGICAL           :: write_xyz    ! writeout xyz files


    LOGICAL                               :: add_parts
    !!! add new particles when needed. Default is .true.
    LOGICAL                               :: del_parts
    !!! delete particles when too many neighbours. Default is .true.
    LOGICAL                               :: remove_large_parts
    !!! delete particles that have a cutoff equal to the maximum cutoff
    !!! This is useful in simulations where not some regions of space
    !!! should be left empty (like in level-sets).  Default is .false.
    LOGICAL                               :: D_needs_gradients
    !!! The monitor function depends on the fields gradient
    REAL(MK)                            :: scale_D        ! resolution scale
    ! maximum scaling distance allowed
    REAL(MK)                            :: maximum_D      
    ! minimum resolution
    REAL(MK)                            :: minimum_D      
    ! stopping criterion for particle adaptation
    REAL(MK)                            :: adaptivity_criterion    
    !counter for the number of gradient descent steps performed
    INTEGER                               :: nb_grad_desc_steps 
    !separation distance below which 2 particles fuse
    REAL(MK)                            :: fuse_radius 

    !some parameters used in various functions
    REAL(MK)                            :: param_a !used in module_funcs
    REAL(MK)                            :: param_d0, param_d1,param_p0 
    REAL(MK)                            :: param_morse ! rho in the Morse
    !potential
    !new particles are generated at distance equal to spawn_radius * D(ip)
    REAL(MK)                            :: spawn_radius

    REAL(MK)                            :: attractive_radius0 
    !distance below which particles attract each other

    REAL(MK)                            :: rcp_over_D

    !====================================================================!
    ! verlet lists
    !====================================================================!
    INTEGER                             :: nneigh_critical
    INTEGER                             :: nneigh_theo
    ! number of neighbours above which something must be wrong
    INTEGER                             :: nneigh_toobig

    !====================================================================!
    ! tolerances in numerical schemes
    !====================================================================!
    REAL(MK)           :: tolerance_moment_conditions
    REAL(MK)           :: tolerance_LSE
    REAL(MK)           :: tolerance_singular_values
    REAL(MK)           :: tolerance_livecheck


END TYPE DTYPE(sop_t_opts)

TYPE DTYPE(sop_t_stats)
    !!! derived type to store information about adaptation steps

    INTEGER                           :: nb_grad_desc_steps
    !!! number of iterations of the gradient descent algorithm
    REAL(MK)                          :: min_sv
    !!! smallest singular value encountered thus far in when computing
    !!! the dc operators (checked only if opts%check_dcops is true)
END TYPE DTYPE(sop_t_stats)

#undef   MK
#undef   _MK
#undef   DTYPE
#undef   CTYPE
