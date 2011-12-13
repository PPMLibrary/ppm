     TYPE DTYPE(sop_t_diff_eq)
         !!! derived type to define a differential operator
         !!! The latter is assumed to be a linear combination
         !!! of multivariate derivatives

         INTEGER                           :: nterms
         !!! number of terms
         REAL(prec),DIMENSION(:),  POINTER :: coeffs
         !!! coefficients in front of each term
         REAL(prec),DIMENSION(:,:),POINTER :: order_d
         !!! degree of the derivative for each term

     END TYPE DTYPE(sop_t_diff_eq)

     TYPE DTYPE(sop_t_opts)
         !!! derived type for optional arguments to the adaptivity routine

         LOGICAL           :: level_set
         REAL(prec)        :: param_nb
         REAL(prec)        :: nb_width
         REAL(prec)        :: nb_width2
         REAL(prec)        :: nb_width_kill
         INTEGER           :: order_approx
         REAL(prec)        :: c
         LOGICAL           :: check_dcops
         !!! order of approximation for the interpolation kernels
         TYPE(DTYPE(sop_t_diff_eq)),POINTER :: diff_eq
         !!! RHS of the differential equation that needs to be computed
         !!! at the interpolation step !TODO (would be convenient...)

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
         REAL(prec)                            :: scale_D        ! resolution scale
         ! maximum scaling distance allowed
         REAL(prec)                            :: maximum_D      
         ! minimum resolution
         REAL(prec)                            :: minimum_D      
         ! stopping criterion for particle adaptation
         REAL(prec)                            :: adaptivity_criterion    
         !counter for the number of gradient descent steps performed
         INTEGER                               :: nb_grad_desc_steps 
         !separation distance below which 2 particles fuse
         REAL(prec)                            :: fuse_radius 

         !some parameters used in various functions
         REAL(prec)                            :: param_a !used in module_funcs
         REAL(prec)                            :: param_d0, param_d1,param_p0 
         REAL(prec)                            :: param_morse ! rho in the Morse
         !potential
         !new particles are generated at distance equal to spawn_radius * D(ip)
         REAL(prec)                            :: spawn_radius

         REAL(prec)                            :: attractive_radius0 
         !distance below which particles attract each other

         REAL(prec)                            :: rcp_over_D

         !====================================================================!
         ! verlet lists
         !====================================================================!
         INTEGER                               :: nneigh_critical
         INTEGER                               :: nneigh_theo
         ! number of neighbours above which something must be wrong
         INTEGER                               :: nneigh_toobig

         !====================================================================!
         ! tolerances in numerical schemes
         !====================================================================!
         REAL(prec)           :: tolerance_moment_conditions
         REAL(prec)           :: tolerance_LSE
         REAL(prec)           :: tolerance_singular_values
         REAL(prec)           :: tolerance_livecheck


     END TYPE DTYPE(sop_t_opts)

     TYPE DTYPE(sop_t_stats)
         !!! derived type to store information about adaptation steps

         INTEGER                           :: nb_grad_desc_steps
         !!! number of iterations of the gradient descent algorithm
         REAL(prec)                        :: min_sv
         !!! smallest singular value encountered thus far in when computing
         !!! the dc operators (checked only if opts%check_dcops is true)
     END TYPE DTYPE(sop_t_stats)

#undef __KIND
#undef DTYPE
#undef prec
#undef _prec
