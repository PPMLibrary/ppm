TYPE DTYPE(ppm_t_part_prop)
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
        PROCEDURE :: print_info => DTYPE(prop_print_info)

END TYPE DTYPE(ppm_t_part_prop)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold neighlists
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_part_prop)
    TYPE(DTYPE(ppm_t_part_prop)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_part_prop)

TYPE DTYPE(ppm_t_neighlist)
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the neighbour list
    INTEGER                                         :: P_id = 0
    !!! Id of the set of particles that this neighbour list refers to
    !!! The default, 0, stands for "self".
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
    LOGICAL                                         :: uptodate = .FALSE.
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin = 0
    !!! smallest number of neighbors on this processor
    INTEGER                                         :: nneighmax = 0
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
    !!! Data structure to describe differential operator
    !!!  It only contains semantic information on the operator
    INTEGER, DIMENSION(:), POINTER                 :: degree =>NULL()
    !!! degree of each term in the linear combination of differential ops 
    INTEGER, DIMENSION(:), POINTER                 :: order =>NULL()
    !!! approximation order of each term 
    REAL(MK), DIMENSION(:), POINTER                :: coeffs =>NULL()
    !!! array where the coefficients in linear combinations of 
    !!! differential ops are stored
    INTEGER                                        :: nterms
    !!! number of terms
    CHARACTER(LEN=ppm_char)                        :: name
    !!! name of the vector-valued property

    CONTAINS
        PRIVATE
        PROCEDURE :: create => DTYPE(desc_create)
        PROCEDURE :: destroy => DTYPE(desc_destroy)

END TYPE DTYPE(ppm_t_opdesc)

TYPE DTYPE(ppm_t_operator)
    !!! Data structure containing all diff operators for a particle set
    !!! 
    TYPE(DTYPE(ppm_t_part_prop)),          POINTER :: ker => NULL()
    !!! where the operators are stored
    TYPE(DTYPE(ppm_t_opdesc)),             POINTER :: desc => NULL()
    !!! small matrices that describe what each operator does
    INTEGER                                        :: P_id = 0
    !!! Id of the set of particles that this operator takes data from.
    !!! The default, 0, stands for "self" (the operator is computed
    !!! on the same set of particles than the one which contains the data).
    INTEGER                                        :: neigh_id = 1
    !!! Id of the neighbour list that should be used
    !!! The default, 1, refers to "self": the list of neighbours within
    !!! the same set of particles. 

    LOGICAL, DIMENSION(ppm_param_length_opsflags)  :: flags
    !!! logical flags
    !!!    ppm_ops_inc_ghosts
    !!!           true if the operator should be computed for ghost 
    !!!           particles too.  Note that the resulting values 
    !!!           will be wrong for the ghost particles
    !!!           that have some neighbours outside the ghost layers. 
    !!!           Default is false.
    !!!    ppm_ops_interp
    !!!          true if the op interpolates data from one set of particles
    !!!    ppm_ops_iscomputed
    !!!          true if the operator has been computed and is uptodate
    !!!    ppm_ops_isdefined
    !!!          true if the operator has been defined 
    !!!    ppm_ops_vector
    !!!          true if each term represents a component (ie the result
    !!!          of the operator should be a vector field, like for gradients)
    !!!          false if the components are added up (like for the divergence)

    CONTAINS
        PRIVATE
        PROCEDURE :: create => DTYPE(op_create)
        PROCEDURE :: destroy => DTYPE(op_destroy)

END TYPE DTYPE(ppm_t_operator)

!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold operators
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_ops)
    TYPE(DTYPE(ppm_t_operator)), POINTER :: t => NULL()
END TYPE DTYPE(ppm_ptr_ops)

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

!----------------------------------------------------------------------
! Container class
!----------------------------------------------------------------------
TYPE, ABSTRACT :: DTYPE(ppm_t_container)
    INTEGER                                         :: nb = 0
    !!! number of different instances
    INTEGER                                         :: iter_id = 0
    !!! index at which the iterator is currently pointing
    INTEGER                                         :: min_id = HUGE(1)
    !!! smallest non-null index in the array
    INTEGER                                         :: max_id = 0
    !!! highest non-null index in the array
    INTEGER                                         :: vec_size = 0
    !!! current allocated size of the array
END TYPE DTYPE(ppm_t_container)

!ABSTRACT INTERFACE
!    FUNCTION DTYPE(begin)(this) RESULT(iterator)
!        IMPORT DTYPE(ppm_t_container)
!        CLASS(DTYPE(ppm_t_container)), INTENT(IN) :: this
!        CLASS(*), POINTER    :: iterator
!    END FUNCTION DTYPE(begin)
!END INTERFACE
!ABSTRACT INTERFACE
!    FUNCTION DTYPE(next)(this) RESULT(iterator)
!        IMPORT DTYPE(ppm_t_container)
!        CLASS(DTYPE(ppm_t_container)), INTENT(IN) :: this
!        CLASS(*), POINTER  :: iterator
!    END FUNCTION DTYPE(next)
!END INTERFACE

!----------------------------------------------------------------------
! Container for properties
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_props)
    TYPE(DTYPE(ppm_ptr_part_prop)),DIMENSION(:),POINTER :: vec => NULL()
    !!! array with pointers to each property for this particle set
    TYPE(DTYPE(ppm_t_part_prop)),              POINTER  :: iterator => NULL()
    CONTAINS
        PROCEDURE       :: begin => DTYPE(begin_prop)
        PROCEDURE       :: next => DTYPE(next_prop)
        PROCEDURE       :: last => DTYPE(last_prop)
        PROCEDURE       :: prev => DTYPE(prev_prop)
        PROCEDURE       :: destroy => DTYPE(prop_container_destroy)
END TYPE DTYPE(ppm_c_props)
!----------------------------------------------------------------------
! Container for operators
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_operators)
    TYPE(DTYPE(ppm_ptr_ops)),DIMENSION(:),POINTER :: vec => NULL()
    !!! array with pointers to each differential operator for this particle set
    TYPE(DTYPE(ppm_t_operator)),          POINTER :: iterator => NULL()
    CONTAINS
        PROCEDURE       :: begin => DTYPE(begin_op)
        PROCEDURE       :: next => DTYPE(next_op)
        PROCEDURE       :: last => DTYPE(last_op)
        PROCEDURE       :: prev => DTYPE(prev_op)
        PROCEDURE       :: destroy => DTYPE(op_container_destroy)
END TYPE DTYPE(ppm_c_operators)
!----------------------------------------------------------------------
! Container for neighbour lists
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_neighlists)
    TYPE(DTYPE(ppm_ptr_t_neighlist)),DIMENSION(:),POINTER :: vec => NULL()
    !!! array with pointers to each neighbor list for this particle set
    TYPE(DTYPE(ppm_t_neighlist)),                 POINTER :: iterator => NULL()
    CONTAINS
        PROCEDURE       :: begin => DTYPE(begin_neigh)
        PROCEDURE       :: next => DTYPE(next_neigh)
        PROCEDURE       :: last => DTYPE(last_neigh)
        PROCEDURE       :: prev => DTYPE(prev_neigh)
        PROCEDURE       :: destroy => DTYPE(neigh_container_destroy)
END TYPE DTYPE(ppm_c_neighlists)

TYPE DTYPE(ppm_t_particles)
    !!! Data structure for a particle set

    CHARACTER(LEN=ppm_char)                         :: name
    !!! name for this particle set
    REAL(MK), DIMENSION(:,:), POINTER               :: xp => NULL()
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
    !!!    ppm_part_neighlists
    !!!          true if neighbour lists within this particle cloud are
    !!!          up-to-date


    ! Topology (links between Particles and topologies)
    INTEGER                                         :: active_topoid


    !ghost layers
    REAL(MK)                                        :: ghostlayer
    !!! ghost layer size used for ghost mapping
    !!! (note: could in principle be different from the cutoff used to 
    !!!  compute the neighbour lists)
    INTEGER                                         :: isymm

    ! Container for Particles properties
    TYPE(DTYPE(ppm_c_props))                        :: props


    ! Container for Neighbor lists
    TYPE(DTYPE(ppm_c_neighlists))                    :: neighs


    ! Container for differential operators
    TYPE(DTYPE(ppm_c_operators))                    :: ops


    ! Load balancing
    REAL(MK), DIMENSION(:), POINTER                 :: pcost=> NULL()
    !!! cost associated to each particle 


    ! stats
    TYPE(DTYPE(particles_stats))                    :: stats
    !!! runtime statistics (e.g. timings, memory)


    ! clocks
    REAL(MK)                                        :: time
    !!! clock for this set of particles
    INTEGER                                         :: itime
    !!! iteration number for this set of particles


    !------------miscellaneous----------

    ! spacings between particles
    REAL(MK)                                      :: h_avg
    !!! (global) average distance between particles
    REAL(MK)                                      :: h_min
    !!! (global) minimum distance between particles
    
    CONTAINS
        PRIVATE
        PROCEDURE     :: create => DTYPE(part_create)
        PROCEDURE     :: destroy => DTYPE(part_destroy)
        PROCEDURE     :: initialize => DTYPE(part_initialize)
        PROCEDURE     :: del_parts => DTYPE(part_del_parts)

        PROCEDURE     :: create_prop => DTYPE(part_prop_create)
        PROCEDURE     :: destroy_prop => DTYPE(part_prop_destroy)

        PROCEDURE     :: create_neighlist => DTYPE(part_neigh_create)
        PROCEDURE     :: destroy_neighlist => DTYPE(part_neigh_destroy)
        PROCEDURE     :: comp_neighlist => DTYPE(part_neighlist)

        PROCEDURE     :: create_op => DTYPE(part_op_create)
        PROCEDURE     :: destroy_op => DTYPE(part_op_destroy)
        PROCEDURE     :: comp_op => DTYPE(part_op_compute)
 !       PROCEDURE     :: updated_positions => DTYPE(part_updated_positions)

        PROCEDURE     :: map_part_push => DTYPE(part_prop_push)
        PROCEDURE     :: map_part_pop => DTYPE(part_prop_pop)
        PROCEDURE     :: map => DTYPE(part_mapping)
        PROCEDURE     :: map_ghosts => DTYPE(part_mapping_ghosts)

        PROCEDURE     :: move => DTYPE(part_move)
        PROCEDURE     :: apply_bc => DTYPE(part_apply_bc)


        PROCEDURE     :: print_info => DTYPE(part_print_info)

        PROCEDURE     :: get_xp => DTYPE(get_xp)
        PROCEDURE     :: set_xp => DTYPE(set_xp)
        PROCEDURE     :: get_dcop => DTYPE(get_dcop)
        PROCEDURE     :: set_dcop => DTYPE(set_dcop)
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

TYPE, EXTENDS(DTYPE(ppm_t_particles)) :: DTYPE(ppm_t_sop)
    !!! an extension of the Particle cloud data structure
    !!! for Self-Organizing Particles

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

!    CONTAINS
!        PRIVATE
!        PROCEDURE     :: create => DTYPE(sop_part_create)


END TYPE DTYPE(ppm_t_sop)

#undef   MK
#undef   _MK
#undef   DTYPE
#undef   CTYPE
