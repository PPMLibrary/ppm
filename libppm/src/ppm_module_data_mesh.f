      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_data_mesh
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_mesh
      !!! This module holds mesh mapping, receive buffers and
      !!! mesh ghost layer mapping buffers.
      !!!
      !!! [NOTE]
      !!! The variables declared in this module should not be accessed by the
      !!! PPM client developer. They are managed interally by the library.

         !----------------------------------------------------------------------
         !  Mesh mapping, send and receive lists
         !----------------------------------------------------------------------

         INTEGER, DIMENSION(:  ), POINTER          :: ppm_mesh_isendfromsub
         !!! list of source subs to send from local processor (local sub number
         !!! on source processor)
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_isendblkstart
         ! start (lower-left corner) of mesh block to be sent in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_isendblksize
         ! size (in grid points) of blocks to be sent
         INTEGER, DIMENSION(:  ), POINTER          :: ppm_mesh_irecvtosub
         ! list of destination subs to recv to on local processors (local sub
         ! number on destination processor)
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_irecvblkstart
         ! start (lower-left corner) of mesh block to be recvd in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_irecvblksize
         ! size (in grid points) of blocks to be recvd



      END MODULE ppm_module_data_mesh
