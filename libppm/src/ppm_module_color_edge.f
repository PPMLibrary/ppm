      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_module_color_edge
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      MODULE ppm_module_color_edge
      !!! Colors edges of a graph using binary heap lists, given number of
      !!! vertices, number of edges, edge list of the graph and optres array
      !!! as input, respectively.

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
          INTEGER, DIMENSION(:), POINTER :: list
          !!! list of vertices that the vertex is connected to
      END TYPE vertex

      TYPE list
      !!! declaration of type: list
          INTEGER, DIMENSION(:), POINTER :: adj_edge
          !!! list of adjacent node of the node
      END TYPE list

      !-------------------------------------------------------------------------
      !  Declaration of arrays
      !-------------------------------------------------------------------------
      TYPE(vertex), ALLOCATABLE, DIMENSION(:) :: node
      !!! Array of nodes
      INTEGER, DIMENSION(:), ALLOCATABLE      :: nelem
      !!! array for number of nodes that are adjacent to each node
      INTEGER, DIMENSION(:), ALLOCATABLE      :: offset
      !!! where to put the next node in the adjacency list of another
      TYPE(list), DIMENSION(:), ALLOCATABLE   :: edges_per_node
      !!! number of edges per node
      TYPE(list), DIMENSION(:), ALLOCATABLE   :: lists
      !!! array of adjacency lists, one for each node
      INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: node_sat
      !!! 2-D array for keeping nodes according to their d-sat values
      !!! and degrees where rows are dsat values and columns are node
      !!! numbers sorted by degree of nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)      :: size_heap
      !!! size of the heap for each row (d-sat value)
      LOGICAL,    ALLOCATABLE, DIMENSION(:)   :: used_color
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
