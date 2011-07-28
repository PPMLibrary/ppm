      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_neighlist_clist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for neighbor
      !                 search routines (cell lists, Verlet lists).
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_neighlist_clist.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:30:02  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschebngraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_neighlist_clist

         !----------------------------------------------------------------------
         !  Define interface to the ppm_neighlist_clist
         !----------------------------------------------------------------------
         INTERFACE ppm_neighlist_clist
            MODULE PROCEDURE ppm_neighlist_clist_d
            MODULE PROCEDURE ppm_neighlist_clist_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_neighlist_clist.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_neighlist_clist.f"
#undef  __KIND

      END MODULE ppm_module_neighlist_clist
