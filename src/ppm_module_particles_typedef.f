MODULE ppm_module_particles_typedef

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __KIND __DOUBLE_PRECISION

USE ppm_module_typedef

IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_double
#endif


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

TYPE pnt_array_1d
    REAL(ppm_kind_double), DIMENSION(:), POINTER :: vec => NULL()
    !!! array where the scalar-value property is stored
    CHARACTER(len=ppm_char)                      :: name
    !!! name of the scalar-valued property
    LOGICAL                                      :: has_ghosts
    !!! true if ghost values are up-to-date
    LOGICAL                                      :: is_mapped
    !!! true if there is a one-to-one mapping with the particles
    LOGICAL                                      :: map_parts
    !!! true if partial mappings are desired for this property (default)
    LOGICAL                                      :: map_ghosts
    !!! true if ghost mappings are desired for this property (default)
END TYPE pnt_array_1d

TYPE pnt_array_2d
    REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: vec =>NULL()
    !!! array where the vector-value property is stored
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the vector-valued property
    INTEGER                                        :: lda
    !!! leading dimension of the array
    LOGICAL                                        :: has_ghosts
    !!! true if ghost values are up-to-date
    LOGICAL                                        :: is_mapped
    !!! true if there is a one-to-one mapping with the particles
    LOGICAL                                      :: map_parts
    !!! true if partial mappings are desired for this property (default)
    LOGICAL                                      :: map_ghosts
    !!! true if ghost mappings are desired for this property (default)
END TYPE pnt_array_2d

TYPE pnt_array_desc
    INTEGER, DIMENSION(:), POINTER                 :: degree =>NULL()
    !!! degree of each term in the linear combination of differential ops 
    INTEGER, DIMENSION(:), POINTER                 :: order =>NULL()
    !!! approximation order of each term 
    REAL(ppm_kind_double), DIMENSION(:), POINTER   :: coeffs =>NULL()
    !!! array where the coefficients in linear combinations of 
    !!! differential ops are stored
    INTEGER                                        :: nterms
    !!! number of terms
    LOGICAL                                        :: vector
    !!! true if each term represents a component (ie. the result
    !!! of the operator should be a vector field, like for the gradient)
    !!! false if the components are added up (like for the divergence)
    CHARACTER(LEN=ppm_char)                        :: name
    !!! name of the vector-valued property
    LOGICAL                                        :: interp
    !!! True if the operator interpolates data from one set of particles
    LOGICAL                                        :: is_computed
    !!! True if the operator has been computed (will be set to False when
    !!! particles move, or e.g. after particle-particle interpolation)
    LOGICAL                                        :: is_defined
    !!! True if the operator has been defined (this usually means that
    !!! the user wants to keep it, unless explicitely stated otherwise)
END TYPE pnt_array_desc

TYPE ppm_t_operator
    TYPE(pnt_array_2d), DIMENSION(:),        POINTER :: ker => NULL()
    !!! array of pointers to where the operators are stored
    TYPE(pnt_array_desc), DIMENSION(:),      POINTER :: desc => NULL()
    !!! array of pointers to small matrices that describe what 
    !!! each operator does
    INTEGER                                          :: nb_ops
    !!! number of operators stored here
    INTEGER                                          :: max_opsid
    !!! largest index for operators
END TYPE ppm_t_operator

TYPE ppm_t_particles

    CHARACTER(LEN=ppm_char)                         :: name
    !!! name for these particles
    REAL(prec), DIMENSION(:,:), POINTER             :: xp => NULL()
    !!! positions of the particles
    LOGICAL                                         :: has_ghosts 
    !!! pseudo-boolean for the positions
    !!! takes value 1 if the ghost values have been computed
    !!! 0 if they have not, and -1 if they do not need to be updated
    INTEGER                                         :: Npart
    !!! Number of real particles on this processor
    INTEGER                                         :: Mpart
    !!! Number of particles (including ghosts) on this processor
    LOGICAL                                         :: areinside
    !!! true if all the particles are inside the comp. domain


    ! Particles properties
    !   integer
    TYPE(pnt_array_1d_i),    DIMENSION(:), POINTER  :: wpi  => NULL()
    !!! array of integer properties
    INTEGER                                         :: nwpi
    !!! number of integer properties
    INTEGER                                         :: max_wpiid
    !!! maximum index for integer properties


    !   scalar
    TYPE(pnt_array_1d),    DIMENSION(:),   POINTER  :: wps  => NULL()
    !!! array of scalar properties
    INTEGER                                         :: nwps
    !!! number of scalar properties
    INTEGER                                         :: max_wpsid
    !!! maximum index for scalar properties


    !   vectors
    TYPE(pnt_array_2d),    DIMENSION(:),   POINTER  :: wpv  => NULL()
    !!! array of vector properties
    INTEGER                                         :: nwpv
    !!! number of vector properties
    INTEGER                                         :: max_wpvid
    !!! maximum index for vector properties


    ! Load balancing
    REAL(prec), DIMENSION(:  ), POINTER             :: pcost=> NULL()
    !!! cost associated to each particle 


    ! Neighbor lists
    REAL(prec)                                      :: cutoff
    !!! cutoff radius
    REAL(prec)                                      :: skin
    !!! skin layer around the particles
    INTEGER                                         :: isymm
    !!! using symmetry
    INTEGER              , DIMENSION(:  ), POINTER  :: nvlist=> NULL()
    !!! Number of neighbors of each particles
    INTEGER              , DIMENSION(:,:), POINTER  :: vlist=> NULL()
    !!! Neighbor lists
    LOGICAL                                         :: neighlists
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin
    !!! smallest number of neighbors
    INTEGER                                         :: nneighmax
    !!! highest number of neighbors

    
    ! xset Neighbor lists
    TYPE(ppm_t_particles),POINTER                   :: Particles_cross=> NULL()
    !!! Pointer to a second set of Particles (not necessary)
    INTEGER              , DIMENSION(:  ), POINTER  :: nvlist_cross=> NULL()
    !!! Number of neighbors of each particles
    INTEGER              , DIMENSION(:,:), POINTER  :: vlist_cross=> NULL()
    !!! Neighbor lists
    LOGICAL                                         :: neighlists_cross
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin_cross
    !!! smallest number of neighbors
    INTEGER                                         :: nneighmax_cross
    !!! highest number of neighbors
    INTEGER                                         :: nn_sq_id
    !!! index of the wps array where nearest-neighbour distances are stored


    ! Links between Particles and topologies
    INTEGER                                         :: active_topoid
    !!! active topology for these particles
    LOGICAL                                         :: ontopology
    !!! true if the particles have been mapped on the active topology


    ! DC-PSE
    INTEGER                                         :: eta_id
    !!! index of the wpv array where the current DC operators are stored
    !!! Will soon be depreciated by this:
    TYPE(ppm_t_operator),  POINTER                  :: ops => NULL()
    !!! structure that contains the discrete operators

    !Global indexing
    INTEGER                                         :: gi_id
    !!! index of the wpi array where the global indices are stored (if needed)

    ! Adaptive particles
    LOGICAL                                         :: adaptive
    !!! true if the particles have their own cutoff radii
    !!! in this case, the cutoff will be stored in wps(rcp_id)%vec
    INTEGER                                         :: rcp_id
    !!! index of the wps array where the cutoff radius is stored
    INTEGER                                         :: D_id
    !!! index of the wps array where D is stored
    INTEGER                                         :: Dtilde_id
    !!! index of the wps array where D_tilde is stored
    INTEGER                                         :: adapt_wpid
    !!! index of the wps array where is stored the property on 
    !!! which adaptation is based 
    !!! default is first 1d property that is not rcp_id (if any)
    !!! otherwise, it is rcp_id
!    INTEGER                                         :: adapt_wpgradid
!    !!! index of the wpv array where is stored the gradient of the property 
!    !!! on which adaptation is based (if needed)
    LOGICAL                                         :: level_set
    !!! true if particles carry a level-set function
    INTEGER                                         :: level_id
    !!! index of the wps array where the level-set is stored
!    INTEGER                                         :: level_old_id
    !!! index of the wps array where the level-set is backed up before adapt

    INTEGER                                         :: level_grad_id
    !!! index of the wps array where the gradient of the level-set is stored
!    INTEGER                                         :: level_grad_old_id
    !!! index of the wps array where the gradient of the level-set 
    !!! is backed up before adapt


    ! Anisotropic particles
    LOGICAL                                         :: anisotropic
    !!! true if the particles have their own cutoff radii
    !!! in this case, the G tensor will be stored in wpv(G_id)%vec
    INTEGER                                         :: G_id
    !!! index where G is stored


    ! spacings between particles
    LOGICAL                                         :: cartesian
    !!! true if particles are located exactly on the nodes of a uniform grid
    REAL(prec)                                      :: h_avg
    !!! (global) average distance between particles
    REAL(prec)                                      :: h_min
    !!! (global) minimum distance between particles

    REAL(prec)                                      :: time
    INTEGER                                         :: itime

END TYPE ppm_t_particles

!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold Particles
!----------------------------------------------------------------------
TYPE ppm_ptr_t_Particles
    TYPE(ppm_t_particles), POINTER  :: t => NULL()
END TYPE ppm_ptr_t_Particles

END MODULE ppm_module_particles_typedef

#undef __KIND