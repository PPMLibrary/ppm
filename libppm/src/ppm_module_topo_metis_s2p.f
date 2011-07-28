      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_topo_metis_s2p
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 decomposition routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_topo_metis_s2p.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 14:12:23  ivos
      !  Renamed from the versions without topo to resolve name conflict
      !  with global ppm_subs2proc data array.
      !
      !  Revision 1.1  2004/07/26 08:55:32  ivos
      !  Renamed.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_topo_metis_s2p

         !----------------------------------------------------------------------
         !  Define interface to METIS-based sub-to-proc assignment
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_metis_s2p
            MODULE PROCEDURE ppm_topo_metis_s2p_s
            MODULE PROCEDURE ppm_topo_metis_s2p_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_topo_metis_s2p.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_topo_metis_s2p.f"
#undef __KIND

      END MODULE ppm_module_topo_metis_s2p
