      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_topo_store
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 callable from the outside. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_topo_store.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:30:10  ivos
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

      MODULE ppm_module_topo_store

         !----------------------------------------------------------------------
         !  Define interfaces to topology store routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_store
            MODULE PROCEDURE ppm_topo_store_s
            MODULE PROCEDURE ppm_topo_store_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_topo_store.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_topo_store.f"
#undef __KIND

      END MODULE ppm_module_topo_store
