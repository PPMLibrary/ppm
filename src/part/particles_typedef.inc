TYPE DTYPE(pnt_array_1d)
    REAL(prec), DIMENSION(:), POINTER              :: vec => NULL()
    !!! array where the scalar-value property is stored
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the scalar-valued property
    LOGICAL                                        :: has_ghosts
    !!! true if ghost values are up-to-date
    LOGICAL                                        :: is_mapped
    !!! true if there is a one-to-one mapping with the particles
    LOGICAL                                        :: map_parts
    !!! true if partial mappings are desired for this property (default)
    !!! (if false, the array for this property is not reallocated when
    !!! particles move to a different processor or when they are
    !!! interpolated from one distribution to another)
    LOGICAL                                        :: map_ghosts
    !!! true if ghost mappings are desired for this property (default)
END TYPE DTYPE(pnt_array_1d)

TYPE DTYPE(pnt_array_2d)
    REAL(prec), DIMENSION(:,:), POINTER            :: vec =>NULL()
    !!! array where the vector-value property is stored
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the vector-valued property
    INTEGER                                        :: lda
    !!! leading dimension of the array
    LOGICAL                                        :: has_ghosts
    !!! true if ghost values are up-to-date
    LOGICAL                                        :: is_mapped
    !!! true if there is a one-to-one mapping with the particles
    LOGICAL                                        :: map_parts
    !!! true if partial mappings are desired for this property (default)
    !!! (if false, the array for this property is not reallocated when
    !!! particles move to a different processor or when they are
    !!! interpolated from one distribution to another)
    LOGICAL                                        :: map_ghosts
    !!! true if ghost mappings are desired for this property (default)
END TYPE DTYPE(pnt_array_2d)

TYPE DTYPE(pnt_array_desc)
    INTEGER, DIMENSION(:), POINTER                 :: degree =>NULL()
    !!! degree of each term in the linear combination of differential ops 
    INTEGER, DIMENSION(:), POINTER                 :: order =>NULL()
    !!! approximation order of each term 
    REAL(prec), DIMENSION(:), POINTER   :: coeffs =>NULL()
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
    !!! true if the operator interpolates data from one set of particles
    LOGICAL                                        :: is_computed
    !!! true if the operator has been computed (will be set to False when
    !!! particles move, or e.g. after particle-particle interpolation)
    LOGICAL                                        :: is_defined
    !!! true if the operator has been defined (this usually means that
    !!! the user wants to keep it, unless explicitely stated otherwise)
    LOGICAL                                        :: with_ghosts
    !!! true if the operator should be computed for ghost particles too. 
    !!! Note that the resulting values will be wrong for the ghost particles
    !!! that have some neighbours outside the ghost layers. Default is false.
END TYPE DTYPE(pnt_array_desc)

TYPE DTYPE(ppm_t_operator)
    TYPE(DTYPE(pnt_array_2d)), DIMENSION(:),  POINTER :: ker => NULL()
    !!! array of pointers to where the operators are stored
    TYPE(DTYPE(pnt_array_desc)), DIMENSION(:),POINTER :: desc => NULL()
    !!! array of pointers to small matrices that describe what 
    !!! each operator does
    INTEGER                                          :: nb_ops
    !!! number of operators stored here
    INTEGER                                          :: max_opsid
    !!! largest index for operators
END TYPE DTYPE(ppm_t_operator)

TYPE DTYPE(particles_stats)
    INTEGER                                          :: nb_nl = 0
    !!! number of neighbour lists constructions
    INTEGER                                          :: nb_inl = 0
    !!! number of inhomogeneous neighbour lists constructions
    INTEGER                                          :: nb_cinl = 0
    !!! number of conventional (old and depreciated) inl constructions
    INTEGER                                          :: nb_xset_inl = 0
    !!! number of xset inhomogeneous neighbour lists constructions
    INTEGER                                          :: nb_xset_nl = 0
    !!! number of xset neighbour lists constructions
    INTEGER                                          :: nb_dc_comp = 0
    !!! number of DC operators computation (matrix inversions)
    INTEGER                                          :: nb_dc_apply = 0
    !!! number of DC operators evaluation (no matrix inversions)
    INTEGER                                          :: nb_kdtree = 0
    !!! number of kdtree constructions
    INTEGER                                          :: nb_global_map = 0
    !!! number of global mappings
    INTEGER                                          :: nb_part_map = 0
    !!! number of partial mappings
    INTEGER                                          :: nb_ghost_get = 0
    !!! number of partial mappings
    INTEGER                                          :: nb_ghost_push = 0
    !!! number of partial mappings
    REAL(prec)                                       :: t_nl = 0._prec
    !!! time spent for neighbour lists constructions
    REAL(prec)                                       :: t_inl = 0._prec
    !!! time spent for inhomogeneous neighbour lists constructions
    REAL(prec)                                       :: t_cinl = 0._prec
    !!! time spent for conventional (old and depreciated) inl constructions
    REAL(prec)                                       :: t_xset_inl = 0._prec
    !!! time spent for xset inhomogeneous neighbour lists constructions
    REAL(prec)                                       :: t_xset_nl = 0._prec
    !!! time spent for xset neighbour lists constructions
    REAL(prec)                                       :: t_dc_comp = 0._prec
    !!! time spent for DC operators computation (matrix inversions)
    REAL(prec)                                       :: t_dc_apply = 0._prec
    !!! time spent for DC operators evaluation (no matrix inversions)
    REAL(prec)                                       :: t_kdtree = 0._prec
    !!! time spent for kdtree constructions
    REAL(prec)                                       :: t_global_map = 0._prec
    !!! time spent for global mappings
    REAL(prec)                                       :: t_part_map = 0._prec
    !!! time spent for partial mappings
    REAL(prec)                                       :: t_ghost_get = 0._prec
    !!! time spent for partial mappings
    REAL(prec)                                       :: t_ghost_push = 0._prec
    !!! time spent for partial mappings

    INTEGER                                          :: nb_ls = 0
    REAL(prec)                                       :: t_ls = 0._prec
    REAL(prec)                                       :: t_add = 0._prec
    REAL(prec)                                       :: t_del = 0._prec
    REAL(prec)                                       :: t_compD = 0._prec

END TYPE DTYPE(particles_stats)

TYPE DTYPE(ppm_t_particles)

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

    INTEGER,                 DIMENSION(:), POINTER  :: wp_t => NULL()
    !!! array containing the type of each property
    INTEGER,                 DIMENSION(:), POINTER  :: idx => NULL()
    !!! array containing the index of each property in the data structure
    !!! corresponding to the type of this property
    


    !   integer
    TYPE(pnt_array_1d_i),    DIMENSION(:), POINTER  :: wpi  => NULL()
    !!! array of integer properties
    INTEGER                                         :: nwpi
    !!! number of integer properties
    INTEGER                                         :: max_wpiid
    !!! maximum index for integer properties


    !   scalar
    TYPE(DTYPE(pnt_array_1d)),    DIMENSION(:),   POINTER :: wps  => NULL()
    !!! array of scalar properties
    INTEGER                                         :: nwps
    !!! number of scalar properties
    INTEGER                                         :: max_wpsid
    !!! maximum index for scalar properties


    !   vectors
    TYPE(DTYPE(pnt_array_2d)),    DIMENSION(:),   POINTER :: wpv  => NULL()
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
    !!! smallest number of neighbors on this processor
    INTEGER                                         :: nneighmax
    !!! highest number of neighbors on this processor
    LOGICAL                                         :: conventionalinl
    !!! true if the inl are computed with the old slow routine (default is
    !!! false)

    
    ! xset Neighbor lists
    TYPE(DTYPE(ppm_t_particles)),POINTER            :: Particles_cross=> NULL()
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
    TYPE(DTYPE(ppm_t_operator)),  POINTER           :: ops => NULL()
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

    ! runtime statistics
    TYPE(DTYPE(particles_stats))                    :: stats

END TYPE DTYPE(ppm_t_particles)

#undef   _prec
#undef   prec
#undef   DEFINE_MK
#undef   DTYPE
#undef   __KIND
