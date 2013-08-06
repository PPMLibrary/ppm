      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                  ppm_module_color_edge
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      MODULE ppm_module_color_edge
      !!! Colors edges of a graph using binary heap lists, given number of
      !!! vertices, number of edges, edge list of the graph and optres array
      !!! as input, respectively.

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Declaration of types
      !-------------------------------------------------------------------------
      TYPE vertex
      !!! declaration of type: vertex
          INTEGER :: degree
          !!! degree of the vertex
          INTEGER :: color
          !!! color of the vertex
          INTEGER :: dsat
          !!! dsat-value of the vertex
          LOGICAL :: iscolored
          !!! TRUE if the vertex is colored
          INTEGER :: loc_heap
          !!! location of vertex in heap list
          INTEGER, DIMENSION(:), POINTER :: list => NULL()
          !!! list of vertices that the vertex is connected to
      END TYPE vertex

      TYPE list
      !!! declaration of type: list
          INTEGER, DIMENSION(:), POINTER :: adj_edge => NULL()
          !!! list of adjacent node of the node
      END TYPE list

      !-------------------------------------------------------------------------
      !  Declaration of arrays
      !-------------------------------------------------------------------------
      TYPE(vertex), DIMENSION(:), POINTER :: node => NULL()
      !!! Array of nodes
      INTEGER, DIMENSION(:), POINTER          :: nelem => NULL()
      !!! array for number of nodes that are adjacent to each node
      INTEGER, DIMENSION(:), POINTER          :: offset => NULL()
      !!! where to put the next node in the adjacency list of another
      TYPE(list), DIMENSION(:), POINTER       :: edges_per_node => NULL()
      !!! number of edges per node
      TYPE(list), DIMENSION(:), POINTER       :: lists => NULL()
      !!! array of adjacency lists, one for each node
      INTEGER, DIMENSION(:,:), POINTER        :: node_sat => NULL()
      !!! 2-D array for keeping nodes according to their d-sat values
      !!! and degrees where rows are dsat values and columns are node
      !!! numbers sorted by degree of nodes
      INTEGER, DIMENSION(:), POINTER          :: size_heap => NULL()
      !!! size of the heap for each row (d-sat value)
      LOGICAL, DIMENSION(:), POINTER          :: used_color => NULL()
      !!! Array to be used to count number of distinct colors

      !-------------------------------------------------------------------------
      !  Declaration of variables
      !-------------------------------------------------------------------------
      INTEGER                                 :: nvertices
      !!! number of vertices in the graph
      INTEGER                                 :: nedges
      !!! number of edges in the graph
      INTEGER                                 :: max_degree
      !!! degree of the graph
      INTEGER                                 :: ncolor
      !!! number of colors to be used
      INTEGER                                 :: alloc_error
      !!! flag for allocation error

      INTERFACE ppm_color_edge
          MODULE PROCEDURE ppm_color_edge
      END INTERFACE

      !-------------------------------------------------------------------------
      !  Privatizing variables and arrays of the module
      !-------------------------------------------------------------------------
      PRIVATE                                 :: node
      PRIVATE                                 :: nelem
      PRIVATE                                 :: offset
      PRIVATE                                 :: edges_per_node
      PRIVATE                                 :: lists
      PRIVATE                                 :: node_sat
      PRIVATE                                 :: size_heap
      PRIVATE                                 :: used_color
      PRIVATE                                 :: nvertices
      PRIVATE                                 :: nedges
      PRIVATE                                 :: max_degree
      PRIVATE                                 :: ncolor
      PRIVATE                                 :: alloc_error

      !-------------------------------------------------------------------------
      !  Location of subroutine to be used for coloring
      !-------------------------------------------------------------------------
      CONTAINS

#include "util/ppm_color_edge.f"

      END MODULE ppm_module_color_edge
