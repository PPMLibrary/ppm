      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_map_connect_send
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the mapping
      !                 routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_connect_send.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:24:33  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 08:24:46  ivos
      !  Check-in after module cleanup. All ppm_map_part_connect* were
      !  renamed to ppm_map_connect*. Also the modules.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_map_connect_send

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:),   POINTER :: id_send,id_inv,csend,crecv
         INTEGER, DIMENSION(:),   POINTER :: sendbuffer,recvbuffer
         INTEGER, DIMENSION(:,:), POINTER :: cd_local,psend

         PRIVATE :: id_send,id_inv,csend,crecv,sendbuffer,recvbuffer,cd_local
         PRIVATE :: psend

         !----------------------------------------------------------------------
         !  Define inteface to ppm_map_connect_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_connect_send
            MODULE PROCEDURE ppm_map_connect_send
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_connect_send.f"

      END MODULE ppm_module_map_connect_send
