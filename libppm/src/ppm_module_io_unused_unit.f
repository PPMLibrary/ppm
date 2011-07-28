      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_io_unused_unit
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the IO routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_io_unused_unit.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:41  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_io_unused_unit

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_unused_unit
         !----------------------------------------------------------------------
         INTERFACE ppm_io_unused_unit
             MODULE PROCEDURE ppm_io_unused_unit
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_io_unused_unit.f"

      END MODULE ppm_module_io_unused_unit
