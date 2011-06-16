      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_data
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Declare global types and variables.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2006/10/10 21:21:16  ivos
      !  Added comment explaining the order of the faces in bcdef.
      !
      !  Revision 1.6  2006/10/10 20:48:43  walther
      !  Added the ppm_ghost_offsets/d.
      !
      !  Revision 1.5  2006/07/04 15:44:43  ivos
      !  Added missing comment.
      !
      !  Revision 1.4  2006/02/03 09:34:02  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.3  2004/07/26 11:48:10  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 10:52:36  hiebers
      !  added ppm_pi
      !
      !  Revision 1.1  2004/07/26 07:28:16  ivos
      !  Initial implementation. Originated from splitting the old ppm
      !  modules.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data

         !----------------------------------------------------------------------
         !  Global TYPEs
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  Header file for global parameters
         !----------------------------------------------------------------------
         INCLUDE 'ppm_param.h'

         !----------------------------------------------------------------------
         !  buffers for communication
         !----------------------------------------------------------------------
         REAL(ppm_kind_double), DIMENSION(  :), POINTER :: &
     &      ppm_sendbufferd, & ! send buffer for particles (double)
     &      ppm_recvbufferd    ! recv buffer for particles (double)

         REAL(ppm_kind_single), DIMENSION(  :), POINTER :: &
     &      ppm_sendbuffers, & ! send buffer for particles (single)
     &      ppm_recvbuffers    ! recv buffer for particles (single)

         INTEGER              , DIMENSION(:,:), POINTER :: &
     &      ppm_ghosthack      ! invert map of ghost for symmetry

         INTEGER              , DIMENSION(  :), POINTER :: &
     &      ppm_psendbuffer, & ! pointer to particles within the send buffer
     &      ppm_precvbuffer    ! pointer to particles within the recv buffer
                               ! both in terms of particle NOT the actual
                               ! position in the buffer

         INTEGER                                        :: &
     &      ppm_nsendbuffer, & ! the size of the send buffer 
     &      ppm_nrecvbuffer    ! the size of the recv buffer
                               ! both in terms of entries in the buffer NOT
                               ! the number of particles 

         INTEGER                                        :: &
     &      ppm_buffer_set     ! the total number of particle fields packed in
                               ! the send buffer, ie. xp, vp is two sets

         ! the original on-processor particle IDs in the order in which
         ! they are in the sendbuffer. Used to push additional particle
         ! data on the buffer in the correct order.
         INTEGER              , DIMENSION(  :), POINTER :: ppm_buffer2part
         INTEGER              , DIMENSION(  :), POINTER :: ppm_buffer_type
         INTEGER              , DIMENSION(  :), POINTER :: ppm_buffer_dim

         !----------------------------------------------------------------------
         !  ghost offset:
         !  ghost particles may have a spatial offset compared to their real
         !  particle - we need to store this offset in terms of the buffer id
         !  in order to be able to push/send/pop ghost coordinates in the 
         !  ghost_get and ghost_put mapping; this is needed when using Verlet
         !  lists: here nothing is remapped and for symmetry we need both to 
         !  put the particle forces etc. and to get the updated particle 
         !  position; for the asymmetric case we do not put the forces back but
         !  we need to get the updated positions of the ghost - NOT getting 
         !  ghosts - because a new mapping is NOT needed but simply reusing
         !  the old mapping push/send/popping the updated coordinates 
         !  The offsets are stored as the pdata in the psendbuffer, i.e.,
         !  ppm_ghost_offset(ibuffer+0) = xp_offset(1,)
         !  ppm_ghost_offset(ibuffer+1) = xp_offset(2,)
         !  ppm_ghost_offset(ibuffer+2) = xp_offset(3,)
         !----------------------------------------------------------------------
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: ppm_ghost_offsets
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_ghost_offsetd

         !----------------------------------------------------------------------
         !  pointers to the subdomains
         !----------------------------------------------------------------------
         ! number of subs on the current processor. index: topoid
         INTEGER              , DIMENSION(:    ), POINTER :: ppm_nsublist
         ! list of subs of the current processor. 1st index: local sub
         ! number. 2nd: topoid
         INTEGER              , DIMENSION(:,:  ), POINTER :: ppm_isublist
         ! total number of subs on all processors. index: topoid
         INTEGER              , DIMENSION(:    ), POINTER :: ppm_nsubs 
         ! extensions of all subs (double and single prec). 1st index:
         ! x,y,(z). 2nd: sub-ID. 3rd: topoid
         REAL(ppm_kind_double), DIMENSION(:,:,:), POINTER :: &
     &      ppm_min_subd,ppm_max_subd
         REAL(ppm_kind_single), DIMENSION(:,:,:), POINTER :: &
     &      ppm_min_subs,ppm_max_subs
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
         ! index 2: the sub id (GLOBAL of ALL the subs!)
         ! index 3: the topoid
         INTEGER              , DIMENSION(:,:,:), POINTER :: ppm_subs_bc

         !----------------------------------------------------------------------
         !  topology
         !----------------------------------------------------------------------
         INTEGER , DIMENSION(:,:), POINTER :: ppm_subs2proc
         ! highest internal topology number (= #of topologies)
         INTEGER                           :: ppm_max_topoid
         ! ID of the current particle topology (in internal numbering)
         INTEGER                           :: ppm_topoid
         ! ID of the current field topology (in internal numbering)
         INTEGER                           :: ppm_field_topoid
         ! user-numbering of the topologies
         INTEGER , DIMENSION(:), POINTER   :: ppm_user_topoid
         ! inverse list: internal numbers indexed by user numbering
         INTEGER , DIMENSION(:), POINTER   :: ppm_internal_topoid
         ! list of neighboring subs of all local subs. 
         !    index 1: neighbor index
         !    index 2: sub id (local index, not global ID!)
         !    index 3: topoid
         INTEGER , DIMENSION(:,:,:), POINTER :: ppm_ineighsubs
         ! number of neighboring subs of all local subs. 
         !    index 1: sub id (local index, not global ID!)
         !    index 2: topoid
         INTEGER , DIMENSION(:,:), POINTER :: ppm_nneighsubs
         ! list of neighboring processors. Index 1: neighbor index, index
         ! 2: topoid
         INTEGER , DIMENSION(:,:), POINTER :: ppm_ineighlist
         ! number of neighboring processors. Index: topoid
         INTEGER , DIMENSION(:  ), POINTER :: ppm_nneighlist
         ! has optimal communication sequence already been determined for a
         ! certain topology (index: topoid)
         LOGICAL , DIMENSION(:  ), POINTER :: ppm_isoptimized
         ! number of communication rounds needed for partial mapping
         ! (index: topoid)
         INTEGER , DIMENSION(:  ), POINTER :: ppm_ncommseq
         ! optimal communication sequence for this processor. 1st index:
         ! communication round, 2nd index: topoid
         INTEGER , DIMENSION(:,:), POINTER :: ppm_icommseq

         ! boundary conditions for the topology
         !   first  index: 1-6 each of the faces: x-, x+, y-, y+, z-, z+
         !   second index: topoid
         INTEGER              , DIMENSION(:,:  ), POINTER :: ppm_bcdef
         ! physical extend of the topology 
         !   first  index: ppm_dim 
         !   second index: topoid
         REAL(ppm_kind_single), DIMENSION(:,:  ), POINTER :: ppm_min_physs, &
        &                                                    ppm_max_physs
         REAL(ppm_kind_double), DIMENSION(:,:  ), POINTER :: ppm_min_physd, &
        &                                                    ppm_max_physd

         !----------------------------------------------------------------------
         !  mapping
         !----------------------------------------------------------------------
         INTEGER                        :: ppm_map_type 
         INTEGER                        :: ppm_nsendlist,ppm_nrecvlist
         INTEGER, DIMENSION(:), POINTER :: ppm_isendlist,ppm_irecvlist

         !----------------------------------------------------------------------
         !  Precision
         !----------------------------------------------------------------------
         INTEGER :: ppm_kind
         INTEGER :: ppm_mpi_kind

         !----------------------------------------------------------------------
         !  Dimensionality
         !----------------------------------------------------------------------
         INTEGER :: ppm_dim

         !----------------------------------------------------------------------
         !  Debugging
         !----------------------------------------------------------------------
         INTEGER :: ppm_debug = 0

         !----------------------------------------------------------------------
         !  Has ppm_init been called?
         !----------------------------------------------------------------------
         LOGICAL :: ppm_initialized = .FALSE.

         !----------------------------------------------------------------------
         !  parallel variables
         !----------------------------------------------------------------------
         INTEGER :: ppm_nproc
         INTEGER :: ppm_rank
         INTEGER :: ppm_comm
         ! relative speeds of the processors (for load balancing)
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_proc_speed

         !----------------------------------------------------------------------
         !  Numerical tolerance. Differences smaller than this are considered
         !  zero.
         !----------------------------------------------------------------------
         REAL(ppm_kind_double) :: ppm_myepsd
         REAL(ppm_kind_single) :: ppm_myepss

         !----------------------------------------------------------------------
         !  Constants (computed in ppm_init)
         !----------------------------------------------------------------------
         REAL(ppm_kind_double) :: ppm_pi_d
         REAL(ppm_kind_single) :: ppm_pi_s

         !----------------------------------------------------------------------
         !  I/O Units
         !----------------------------------------------------------------------
         INTEGER               :: ppm_stdout = 6
         INTEGER               :: ppm_stderr = 0
         INTEGER               :: ppm_logfile = -1
         
      END MODULE ppm_module_data
