      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_map_part_send
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
      !  $Log: ppm_module_map_part_send.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2004/12/02 10:02:32  ivos
      !  bugfix: ppm_kind_* are now declared PRIVATE to avoid name clashes
      !  with user programs that include ppm_param.h.
      !
      !  Revision 1.2  2004/11/11 15:26:19  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 07:29:56  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map_part_send

         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: sends,recvs
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: sendd,recvd
         INTEGER, DIMENSION(:), POINTER   :: nsend,nrecv,psend,precv
         INTEGER, DIMENSION(:,:), POINTER :: pp,qq

         PRIVATE :: sends,recvs,sendd,recvd,nsend,nrecv,psend,precv,pp,qq

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_send
            MODULE PROCEDURE ppm_map_part_send
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_part_send.f"

      END MODULE ppm_module_map_part_send
