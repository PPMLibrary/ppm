TYPE DTYPE(part_prop)
    !!! Data structure for particle properties
    INTEGER                                        :: data_type
    !!! type of the property
    !!! One of:
    !!!     ppm_param_...
    !!! 
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the property
    LOGICAL, DIMENSION(ppm_param_length_pptflags)  :: flags
    !!! logical flags
    !!!    ppm_ppt_ghosts
    !!!          true if ghost values are up-to-date
    !!!    ppm_ppt_partial
    !!!          true if there is a one-to-one mapping with the particles
    !!!    ppm_ppt_reqput
    !!!    ppm_ppt_map_parts
    !!!          true if partial mappings are desired for this property (default)
    !!!          (if false, the array for this property is not reallocated when
    !!!           particles move to a different processor or when they are
    !!!           interpolated from one distribution to another)
    !!!    ppm_ppt_map_ghosts
    !!!          true if ghost mappings are desired for this property (default)
    INTEGER                                        :: lda
    !!! leading dimension of the data array
    !!!
    !!! pointers to arrays where the scalar-value properties are stored
    INTEGER, DIMENSION(:), POINTER                 :: data_1d_i =>NULL()
    !!! if the data is 1d integer
    INTEGER, DIMENSION(:,:), POINTER               :: data_2d_i =>NULL()
    !!! if the data is 2d integer
    INTEGER(ppm_kind_int64),DIMENSION(:), POINTER  :: data_1d_li =>NULL()
    !!! if the data is 1d long integer
    INTEGER(ppm_kind_int64),DIMENSION(:,:),POINTER :: data_2d_li =>NULL()
    !!! if the data is 2d long integer
    REAL(MK), DIMENSION(:), POINTER                :: data_1d_r => NULL()
    !!! if the data is 1d real
    REAL(MK), DIMENSION(:,:), POINTER              :: data_2d_r =>NULL()
    !!! if the data is 2d real
    COMPLEX(MK), DIMENSION(:), POINTER             :: data_1d_c => NULL()
    !!! if the data is 1d complex
    COMPLEX(MK), DIMENSION(:,:), POINTER           :: data_2d_c =>NULL()
    !!! if the data is 2d complex
    LOGICAL, DIMENSION(:), POINTER                 :: data_1d_l => NULL()
    !!! if the data is 1d logical
    LOGICAL, DIMENSION(:,:), POINTER               :: data_2d_l =>NULL()
    !!! if the data is 2d logical

    CONTAINS
        PRIVATE
        PROCEDURE :: create => DTYPE(prop_create)
        PROCEDURE :: destroy => DTYPE(prop_destroy)

END TYPE DTYPE(part_prop)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold neighlists
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_part_prop)
    TYPE(DTYPE(part_prop)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_part_prop)

TYPE DTYPE(ppm_t_neighlist)
    INTEGER                                         :: P_id
    !!! Id of the set of particles that this neighbour list refers to
    REAL(MK)                                        :: cutoff
    !!! cutoff radius
    REAL(MK)                                        :: skin
    !!! skin layer around the particles
    INTEGER                                         :: isymm
    !!! using symmetry
    INTEGER              , DIMENSION(:  ), POINTER  :: nvlist=> NULL()
    !!! Number of neighbors of each particles
    INTEGER              , DIMENSION(:,:), POINTER  :: vlist=> NULL()
    !!! Neighbor lists
    LOGICAL                                         :: uptodate
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin
    !!! smallest number of neighbors on this processor
    INTEGER                                         :: nneighmax
    !!! highest number of neighbors on this processor

    CONTAINS
        PRIVATE
        PROCEDURE :: destroy => DTYPE(neigh_destroy)

END TYPE DTYPE(ppm_t_neighlist)

!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold neighlists
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_t_neighlist)
    TYPE(DTYPE(ppm_t_neighlist)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_t_neighlist)

TYPE DTYPE(ppm_t_opdesc)
    !!! Data structure for differential operator
    !!!
    INTEGER, DIMENSION(:), POINTER                 :: degree =>NULL()
    !!! degree of each term in the linear combination of differential ops 
    INTEGER, DIMENSION(:), POINTER                 :: order =>NULL()
    !!! approximation order of each term 
    REAL(MK), DIMENSION(:), POINTER                :: coeffs =>NULL()
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

    CONTAINS
        PRIVATE
        PROCEDURE :: destroy => DTYPE(desc_destroy)

END TYPE DTYPE(ppm_t_opdesc)

!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold op descriptions
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_opdesc)
    TYPE(DTYPE(ppm_t_opdesc)), POINTER :: t => NULL()
END TYPE DTYPE(ppm_ptr_opdesc)


TYPE DTYPE(ppm_t_operator)
    !!! Data structure containing all diff operators for a particle set
    !!! 
    TYPE(DTYPE(ppm_ptr_part_prop)),DIMENSION(:),POINTER :: ker => NULL()
    !!! where the operators are stored
    TYPE(DTYPE(ppm_ptr_opdesc)), DIMENSION(:),  POINTER :: desc => NULL()
    !!! small matrices that describe what each operator does
    INTEGER                                          :: nb_ops
    !!! number of operators stored here
    INTEGER                                          :: max_opsid
    !!! largest index for operators

    CONTAINS
        PRIVATE
        PROCEDURE :: destroy => DTYPE(op_destroy)

END TYPE DTYPE(ppm_t_operator)

TYPE DTYPE(particles_stats)
    !!! Data structure containing statistics about a particle set
    !!! 
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
    REAL(MK)                                       :: t_nl = 0._MK
    !!! time spent for neighbour lists constructions
    REAL(MK)                                       :: t_inl = 0._MK
    !!! time spent for inhomogeneous neighbour lists constructions
    REAL(MK)                                       :: t_cinl = 0._MK
    !!! time spent for conventional (old and depreciated) inl constructions
    REAL(MK)                                       :: t_xset_inl = 0._MK
    !!! time spent for xset inhomogeneous neighbour lists constructions
    REAL(MK)                                       :: t_xset_nl = 0._MK
    !!! time spent for xset neighbour lists constructions
    REAL(MK)                                       :: t_dc_comp = 0._MK
    !!! time spent for DC operators computation (matrix inversions)
    REAL(MK)                                       :: t_dc_apply = 0._MK
    !!! time spent for DC operators evaluation (no matrix inversions)
    REAL(MK)                                       :: t_kdtree = 0._MK
    !!! time spent for kdtree constructions
    REAL(MK)                                       :: t_global_map = 0._MK
    !!! time spent for global mappings
    REAL(MK)                                       :: t_part_map = 0._MK
    !!! time spent for partial mappings
    REAL(MK)                                       :: t_ghost_get = 0._MK
    !!! time spent for partial mappings
    REAL(MK)                                       :: t_ghost_push = 0._MK
    !!! time spent for partial mappings

    INTEGER                                          :: nb_ls = 0
    REAL(MK)                                       :: t_ls = 0._MK
    REAL(MK)                                       :: t_add = 0._MK
    REAL(MK)                                       :: t_del = 0._MK
    REAL(MK)                                       :: t_compD = 0._MK

END TYPE DTYPE(particles_stats)

TYPE DTYPE(ppm_t_particles)
    !!! Data structure for a particle set

    CHARACTER(LEN=ppm_char)                         :: name
    !!! name for this particle set
    REAL(MK), DIMENSION(:,:), POINTER             :: xp => NULL()
    !!! positions of the particles
    INTEGER                                         :: Npart
    !!! Number of real particles on this processor
    INTEGER                                         :: Mpart
    !!! Number of particles (including ghosts) on this processor
    LOGICAL, DIMENSION(ppm_param_length_partflags)  :: flags
    !!! logical flags
    !!!    ppm_part_ghosts
    !!!          true if ghost values are up-to-date
    !!!    ppm_part_partial
    !!!          true if the particles have been mapped on the active topology
    !!!    ppm_part_reqput
    !!!    ppm_part_areinside
    !!!          true if all the particles are inside the comp. domain
    !!!    ppm_part_cartesian
    !!!          true if the particles form a Cartesian grid


    ! Topology (links between Particles and topologies)
    INTEGER                                         :: active_topoid

    ! Particles properties
    TYPE(DTYPE(ppm_ptr_part_prop)),DIMENSION(:),POINTER :: props => NULL()
    !!! array with pointers to each property for this particle set
    INTEGER                                         :: nwp
    !!! number of properties
    INTEGER                                         :: max_wpid
    !!! highest non-null index in the props array
    INTEGER                                         :: swp
    !!! current allocated size of the props array


    ! Load balancing
    REAL(MK), DIMENSION(:), POINTER                 :: pcost=> NULL()
    !!! cost associated to each particle 


    ! Neighbor lists
    TYPE(DTYPE(ppm_ptr_t_neighlist)),DIMENSION(:),POINTER :: neighs => NULL()
    !!! array with pointers to each neighbor list for this particle set
    INTEGER                                         :: nnl
    !!! number of different neighbor lists
    INTEGER                                         :: max_nlid
    !!! highest non-null index in the neighs array
    INTEGER                                         :: snl
    !!! current allocated size of the neighs array


    ! DC-PSE
    TYPE(DTYPE(ppm_t_operator)),  POINTER           :: ops => NULL()
    !!! structure that contains the discrete operators


    ! stats
    TYPE(DTYPE(particles_stats))                    :: stats
    !!! runtime statistics (e.g. timings, memory)


    ! clocks
    REAL(MK)                                        :: time
    !!! clock for this set of particles
    INTEGER                                         :: itime
    !!! iteration number for this set of particles


    !------------miscellaneous----------
    INTEGER                                         :: nn_sq_id
    !!! index of the wps array where nearest-neighbour distances are stored

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
    REAL(MK)                                      :: h_avg
    !!! (global) average distance between particles
    REAL(MK)                                      :: h_min
    !!! (global) minimum distance between particles
    
    CONTAINS
        PRIVATE
        PROCEDURE     :: create => DTYPE(part_create)
        PROCEDURE     :: destroy => DTYPE(part_destroy)
        PROCEDURE     :: initialize => DTYPE(part_initialize)
        PROCEDURE     :: get_xp => DTYPE(get_xp)
        PROCEDURE     :: set_xp => DTYPE(set_xp)
        PROCEDURE     :: DTYPE(data_1d_i_get)
        PROCEDURE     :: DTYPE(data_2d_i_get)
        PROCEDURE     :: DTYPE(data_1d_li_get)
        PROCEDURE     :: DTYPE(data_2d_li_get)
        PROCEDURE     :: DTYPE(data_1d_r_get)
        PROCEDURE     :: DTYPE(data_2d_r_get)
        PROCEDURE     :: DTYPE(data_1d_c_get)
        PROCEDURE     :: DTYPE(data_2d_c_get)
        PROCEDURE     :: DTYPE(data_1d_l_get)
        PROCEDURE     :: DTYPE(data_2d_l_get)

        PROCEDURE     :: DTYPE(data_1d_i_set)
        PROCEDURE     :: DTYPE(data_2d_i_set)
        PROCEDURE     :: DTYPE(data_1d_li_set)
        PROCEDURE     :: DTYPE(data_2d_li_set)
        PROCEDURE     :: DTYPE(data_1d_r_set)
        PROCEDURE     :: DTYPE(data_2d_r_set)
        PROCEDURE     :: DTYPE(data_1d_c_set)
        PROCEDURE     :: DTYPE(data_2d_c_set)
        PROCEDURE     :: DTYPE(data_1d_l_set)
        PROCEDURE     :: DTYPE(data_2d_l_set)

        PROCEDURE     :: create_prop => DTYPE(part_prop_create)
        PROCEDURE     :: destroy_prop => DTYPE(part_prop_destroy)

        GENERIC, PUBLIC :: get =>  &
            DTYPE(data_1d_i_get),&
            DTYPE(data_2d_i_get),&
            DTYPE(data_1d_li_get),&
            DTYPE(data_2d_li_get),&
            DTYPE(data_1d_r_get),&
            DTYPE(data_2d_r_get),&
            DTYPE(data_1d_c_get),&
            DTYPE(data_2d_c_get),&
            DTYPE(data_1d_l_get),&
            DTYPE(data_2d_l_get)
        GENERIC, PUBLIC :: set =>  &
            DTYPE(data_1d_i_set),&
            DTYPE(data_2d_i_set),&
            DTYPE(data_1d_li_set),&
            DTYPE(data_2d_li_set),&
            DTYPE(data_1d_r_set),&
            DTYPE(data_2d_r_set),&
            DTYPE(data_1d_c_set),&
            DTYPE(data_2d_c_set),&
            DTYPE(data_1d_l_set),&
            DTYPE(data_2d_l_set)

END TYPE DTYPE(ppm_t_particles)

#undef   MK
#undef   _MK
#undef   DTYPE
#undef   CTYPE
