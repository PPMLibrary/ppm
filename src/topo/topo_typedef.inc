         TYPE DTYPE(ppm_t_topo)
         !!! The topology type

            INTEGER                                      :: ID
            !!! ID of this topology
            !!!
            !!! It is the same as its index in the ppm_topo array
            LOGICAL                                      :: isdefined = .FALSE.
            !!! flag to tell if this topology is defined/in use
            INTEGER                                      :: prec
            !!! numerical precision (ppm_kind) for this topology

            REAL(MK), DIMENSION(:), POINTER                :: min_phys => NULL()
            !!! minimum of physical extend of the computational domain
            !!! Note: first index is ppm_dim
            REAL(MK), DIMENSION(:), POINTER                :: max_phys => NULL()
            !!! maximum of physical extend of the computational domain
            !!! Note: first index is ppm_dim

            INTEGER              , DIMENSION(:  ), POINTER :: bcdef => NULL()
            !!! boundary conditions for the topology
            !!! Note: first index is 1-6 (each of the faces)

            INTEGER                                        :: nsubs
            !!! total number of subs on all processors.

            REAL(MK), DIMENSION(:,:), POINTER :: min_sub => NULL()
            !!! mimimum of extension of subs
            !!! Note: 1st index: x,y,(z), 2nd: subID
            REAL(MK), DIMENSION(:,:), POINTER :: max_sub => NULL()
            !!! maximum of extension of subs
            !!! Note: 1st index: x,y,(z), 2nd: subID

            REAL(MK), DIMENSION(:  ),POINTER :: sub_cost => NULL()
            !!! estimated cost associated with subdomains. Index: sub-ID.

            INTEGER              , DIMENSION(:  ),POINTER :: sub2proc => NULL()
            !!! subdomain to processor assignment. index: subID (global)

            INTEGER                                       :: nsublist
            !!! number of subs on the current processor.

            INTEGER              , DIMENSION(:  ),POINTER :: isublist => NULL()
            !!! list of subs of the current processor. 1st index: local sub
            !!! number.

            INTEGER              , DIMENSION(:,:),POINTER :: subs_bc => NULL()
            !!! boundary conditions on a sub:
            !!!
            !!! - west  : 1
            !!! - east  : 2
            !!! - south : 3
            !!! - north : 4
            !!! - bottom: 5
            !!! - top   : 6
            !!!
            !!! index 1: the index of the 4 or 6 faces in 2 and 3 D
            !!! index 2: the global sub id
            !!!
            !!! states:
            !!!
            !!! - value: 0 the face is internal
            !!! - value: 1 otherwise

            INTEGER            , DIMENSION(:,:), POINTER :: ineighsubs => NULL()
            !!! list of neighboring subs of all local subs.
            !!! - index 1: neighbor index
            !!! - index 2: sub id (local index, not global ID!)

            INTEGER            , DIMENSION(:  ), POINTER :: nneighsubs => NULL()
            !!! number of neighboring subs of all local subs.
            !!!
            !!! index 1: sub id (local index, not global ID!)
            INTEGER            , DIMENSION(:  ), POINTER :: ineighproc => NULL()
            !!! list of neighboring processors. Index 1: neighbor index
            INTEGER                                        :: nneighproc
            !!! number of neighboring processors.
            LOGICAL                                        :: isoptimized
            !!! has optimal communication sequence already been determined for
            !!! this topology?
            INTEGER                                        :: ncommseq
            !!! number of communication rounds needed for partial mapping
            INTEGER            , DIMENSION(:  ), POINTER :: icommseq => NULL()
            !!! optimal communication sequence for this processor. 1st index:
            !!! communication round
            INTEGER                                        :: max_meshid
            !!! Number of meshes defined on this topology

            TYPE(DTYPE(ppm_t_equi_mesh)), DIMENSION(:  ), POINTER :: mesh => NULL()
            !!! List of meshes defined on this topology. Index: meshid

            REAL(MK)                                       :: ghostsize
            !!! max ghostsize width used when creating this topology
            !!! using a larger cutoff when calling routines such as ghost_get
            !!! should ideally raise a warning
         END TYPE DTYPE(ppm_t_topo)

         !----------------------------------------------------------------------
         ! Wrapper type to be able to have a pointer array to hold topologies
         !----------------------------------------------------------------------
         TYPE DTYPE(ppm_ptr_t_topo)
             TYPE(DTYPE(ppm_t_topo)), POINTER  :: t => NULL()
         END TYPE DTYPE(ppm_ptr_t_topo)
         
         

         !----------------------------------------------------------------------
         !  Operator interfaces
         !----------------------------------------------------------------------
!         INTERFACE OPERATOR (=)
!             SUBROUTINE ppm_topo_copy(t1,t2,info)
!                 IMPLICIT NONE
!                 TYPE(ppm_t_topo), INTENT(IN   ) :: t1
!                 TYPE(ppm_t_topo), INTENT(  OUT) :: t2
!                 INTEGER         , INTENT(  OUT) :: info
!             END SUBROUTINE
!         END INTERFACE
!         INTERFACE OPERATOR (=)
!             SUBROUTINE ppm_mesh_copy(m1,m2,info)
!                 IMPLICIT NONE
!                 TYPE(ppm_t_equi_mesh), INTENT(IN   ) :: m1
!                 TYPE(ppm_t_equi_mesh), INTENT(  OUT) :: m2
!                 INTEGER              , INTENT(  OUT) :: info
!             END SUBROUTINE
!         END INTERFACE



#undef DTYPE
#undef MK
