      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_mktopo
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 that create topologies.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mktopo.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/09/24 15:02:04  ivos
      !  Added ppm_topo_mkgeom and ppm_topo_mktree.
      !
      !  Revision 1.1  2004/07/26 07:30:01  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
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

      MODULE ppm_module_mktopo

         !----------------------------------------------------------------------
         !  Define interfaces to the main topology routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_mktopo
            MODULE PROCEDURE ppm_topo_mkpart_s
            MODULE PROCEDURE ppm_topo_mkpart_d
            MODULE PROCEDURE ppm_topo_mkfield_s
            MODULE PROCEDURE ppm_topo_mkfield_d
            MODULE PROCEDURE ppm_topo_mkgeom_s
            MODULE PROCEDURE ppm_topo_mkgeom_d
            MODULE PROCEDURE ppm_topo_mktree_s
            MODULE PROCEDURE ppm_topo_mktree_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_topo_mkpart.f"
#include "ppm_topo_mkfield.f"
#include "ppm_topo_mkgeom.f"
#include "ppm_topo_mktree.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_topo_mkpart.f"
#include "ppm_topo_mkfield.f"
#include "ppm_topo_mkgeom.f"
#include "ppm_topo_mktree.f"
#undef __KIND

      END MODULE ppm_module_mktopo
