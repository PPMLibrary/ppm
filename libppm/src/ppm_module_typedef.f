      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_typedef
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Declare ppm data types accessible to user.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_typedef.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:07  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_typedef

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data
         !----------------------------------------------------------------------
         !  TYPE for equispaced cartesian meshes on subs
         !----------------------------------------------------------------------
         TYPE ppm_t_equi_mesh
             ! The number of mesh NODES (not cells) in each direction in
             ! each sub
             INTEGER, DIMENSION(:,:), POINTER  :: nnodes
             ! Starting indices of the mesh of this sub in the global mesh
             INTEGER, DIMENSION(:,:), POINTER  :: istart
             ! global number of mesh points in computational domain
             INTEGER, DIMENSION(:  ), POINTER  :: Nm
         END TYPE

         !----------------------------------------------------------------------
         !  Topology TYPE
         !----------------------------------------------------------------------
         TYPE ppm_t_topo
            ! ID (handle) of this topology
            INTEGER                                        :: ID
            ! flag to tell if this topology is defined/in use
            LOGICAL                                        :: isdefined
            ! numerical precision (ppm_kind) for this topology
            INTEGER                                        :: prec
            ! physical extend of the computational domain
            !   first  index: ppm_dim 
            REAL(ppm_kind_single), DIMENSION(:), POINTER :: min_physs,max_physs
            REAL(ppm_kind_double), DIMENSION(:), POINTER :: min_physd,max_physd
            ! boundary conditions for the topology
            !   first  index: 1-6 each of the faces
            INTEGER              , DIMENSION(:  ), POINTER :: bcdef
            ! total number of subs on all processors.
            INTEGER                                        :: nsubs
            ! extensions of all subs (double and single prec). 1st index:
            ! x,y,(z). 2nd: sub-ID.
            REAL(ppm_kind_single), DIMENSION(:,:), POINTER :: min_subs,max_subs
            REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: min_subd,max_subd
            ! estimated cost associated with subdomains. Index: sub-ID.
            REAL(ppm_kind_single), DIMENSION(:  ), POINTER :: sub_costs
            REAL(ppm_kind_double), DIMENSION(:  ), POINTER :: sub_costd
            ! subdomain to processor assignment. index: subID (global)
            INTEGER              , DIMENSION(:  ), POINTER :: subs2proc
            ! number of subs on the current processor.
            INTEGER                                        :: nsublist
            ! list of subs of the current processor. 1st index: local sub
            ! number.
            INTEGER              , DIMENSION(:  ), POINTER :: isublist
            ! boundary conditions on a sub:
            !    west  : 1
            !    east  : 2
            !    south : 3
            !    north : 4
            !    bottom: 5
            !    top   : 6
            ! index 1: the state of the 4 or 6 faces in 2 and 3 D
            !    value: 0 the face is internal
            !    value: 1 otherwise
            ! index 2: the local sub id (ONLY the local subs !)
            INTEGER              , DIMENSION(:,:), POINTER :: subs_bc
            ! list of neighboring subs of all local subs. 
            !    index 1: neighbor index
            !    index 2: sub id (local index, not global ID!)
            INTEGER              , DIMENSION(:,:), POINTER :: ineighsubs
            ! number of neighboring subs of all local subs. 
            !    index 1: sub id (local index, not global ID!)
            INTEGER              , DIMENSION(:  ), POINTER :: nneighsubs
            ! list of neighboring processors. Index 1: neighbor index
            INTEGER              , DIMENSION(:  ), POINTER :: ineighproc
            ! number of neighboring processors. 
            INTEGER                                        :: nneighproc
            ! has optimal communication sequence already been determined for
            ! this topology?
            LOGICAL                                        :: isoptimized
            ! number of communication rounds needed for partial mapping
            INTEGER                                        :: ncommseq
            ! optimal communication sequence for this processor. 1st index:
            ! communication round
            INTEGER              , DIMENSION(:  ), POINTER :: icommseq
            ! Number of meshes defined on this topology
            INTEGER                                        :: max_meshid
            ! List of meshes defined on this topology. Index: meshid
            TYPE(ppm_t_equi_mesh), DIMENSION(:  ), POINTER :: mesh
         END TYPE ppm_t_topo

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

      END MODULE ppm_module_typedef
