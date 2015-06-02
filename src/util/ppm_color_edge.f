      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_color_edge
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
      SUBROUTINE ppm_color_edge(numV,edge_array,coloring,info)
      !!! Given the edge array as input and coloring array to be
      !!! modified, colors edges and updates coloring array such that
      !!! coloring array looks like (p1,p2,c1, ..., pX,pY,cZ)
      !!!
      !!! [NOTE]
      !!! This subroutine was introduced in PPM library to replace coloring
      !!! algorithm of Vizing which was in C++

      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,               INTENT(IN   ) :: numV
      INTEGER, DIMENSION(:), INTENT(INOUT) :: edge_array
      INTEGER, DIMENSION(:), INTENT(INOUT) :: coloring
      INTEGER,               INTENT(  OUT) :: info

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: i,iopt
      INTEGER :: idx

      CHARACTER(LEN=ppm_char) :: caller="ppm_color_edge"

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      nvertices = numV
      !-------------------------------------------------------------------------
      !  Using the edge array provided, construct adjacency lists and assign
      !  them to corresponding nodes
      !-------------------------------------------------------------------------
      CALL create_adjacency_lists(edge_array,info)
      or_fail("create_adjacency_lists")
      !-------------------------------------------------------------------------
      !  Allocate memory for binary heap lists and initialize them
      !-------------------------------------------------------------------------
      CALL initialize_binary_heap(info)
      or_fail("initialize_binary_heap")
      !-------------------------------------------------------------------------
      !  Initialize nodes by setting their variables to prior values
      !-------------------------------------------------------------------------
      CALL initialize_nodes(info)
      or_fail("initialize_nodes")

      !-------------------------------------------------------------------------
      !  Insert all nodes in binary heap lists. They will all go to 0th row of
      !  binary heap lists as all of them has dsat-value 0 in the beginning
      !-------------------------------------------------------------------------
      DO i = 1,nedges
         CALL insert_node(i,info)
         or_fail("insert_node")
      ENDDO

      !-------------------------------------------------------------------------
      !  All edges are colored sequentially, by getting the node to color,
      !  coloring the node and then, updating its neighbors so that they will be
      !  shifted to correct row of binary heap lists if necessary
      !-------------------------------------------------------------------------
      DO i = 1, nedges
         idx = next_node_to_color()
         CALL color_edge(idx,info)
         or_fail("color_edge")

         CALL update_neighbors(idx,info)
         or_fail("update_neighbors")
      ENDDO

      !-------------------------------------------------------------------------
      !  coloring array is modified such that it is of the form (p1,p2,c) ...
      !-------------------------------------------------------------------------
      DO i = 1, nedges
         coloring(3*i-2) = edge_array(2*i-1)
         coloring(3*i-1) = edge_array(2*i)
         coloring(3*i)   = node(i)%color
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate everything that was used
      !-------------------------------------------------------------------------
      DEALLOCATE(used_color,STAT=info)
      or_fail_dealloc("used_color")

      iopt=ppm_param_dealloc
      DO i = 1, nedges
         CALL ppm_alloc(lists(i)%adj_edge,ldc,iopt,info)
         or_fail_dealloc("lists(i)%adj_edge")

         NULLIFY(node(i)%list)
      ENDDO

      DEALLOCATE(node,lists,node_sat,size_heap,STAT=info)
      or_fail_dealloc("node,lists,node_sat & size_heap")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)

      CONTAINS
          !---------------------------------------------------------------------
          ! Subroutine that calls everything in order and
          ! DEALLOCATEs intermediate lists
          !---------------------------------------------------------------------
          SUBROUTINE create_adjacency_lists(input_array,info)

          IMPLICIT NONE

          INTEGER, DIMENSION(:), INTENT(INOUT) :: input_array
          INTEGER,               INTENT(  OUT) :: info

          INTEGER :: i,iopt

          CHARACTER(LEN=ppm_char) :: caller="create_adjacency_lists"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          CALL order_vertices(input_array,info)
          or_fail("order_vertices")

          CALL allocate_processor_lists(input_array,info)
          or_fail("allocate_processor_lists")

          CALL fill_processor_lists(input_array,info)
          or_fail("fill_processor_lists")

          CALL allocate_edge_lists(input_array,info)
          or_fail("allocate_edge_lists")

          CALL fill_edge_lists(input_array,info)
          or_fail("fill_edge_lists")

          CALL assign_edge_lists(info)
          or_fail("assign_edge_lists")

          CALL get_maximum_degree(info)
          or_fail("get_maximum_degree")

          DEALLOCATE(nelem,offset,STAT=info)
          or_fail_dealloc("nelem & offset")

          iopt=ppm_param_dealloc
          DO i = 1,nvertices
             CALL ppm_alloc(edges_per_node(i)%adj_edge,ldc,iopt,info)
             or_fail_dealloc("edges_per_node(i)%adj_edge")
          ENDDO

          DEALLOCATE(edges_per_node,STAT=info)
          or_fail_dealloc("edges_per_node")
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          9999 CONTINUE
          CALL substop(caller,t0,info)
          END SUBROUTINE create_adjacency_lists

          !---------------------------------------------------------------------
          !  Given an array of edges, sorts pairs such that first value is
          !  smaller than the second, to guarantee e1<e2 at all times
          !---------------------------------------------------------------------
          SUBROUTINE order_vertices(input_array,info)

          IMPLICIT NONE

          INTEGER, DIMENSION(:), INTENT(INOUT) :: input_array
          INTEGER,               INTENT(  OUT) :: info

          INTEGER :: i
          INTEGER :: dummy

          start_subroutine("order_vertices")

          nedges = size(input_array)/2
          DO i = 1,nedges
             IF (input_array(2*i-1).GT.input_array(2*i)) THEN !swap vertices
                dummy = input_array(2*i-1)
                input_array(2*i-1) = input_array(2*i)
                input_array(2*i) = dummy
             ENDIF
          ENDDO

          end_subroutine()
          END SUBROUTINE order_vertices

          !---------------------------------------------------------------------
          !  Allocates adjacency lists for every processor, a list for
          !  each processor
          !---------------------------------------------------------------------
          SUBROUTINE allocate_processor_lists(input_array,info)

          IMPLICIT NONE

          INTEGER, DIMENSION(:), INTENT(IN   ) :: input_array
          INTEGER,               INTENT(  OUT) :: info

          INTEGER :: i,iopt

          start_subroutine("allocate_processor_lists")

          ALLOCATE(nelem(1:nvertices),STAT=info)
          or_fail_alloc("Could not allocate nelem")

          ALLOCATE(edges_per_node(1:nvertices),STAT=info)
          or_fail_alloc("Could not allocate edges_per_node")

          nelem = 0
          !count number of processors that are connected
          DO i = 1,nedges
             nelem(input_array(2*i-1)) = nelem(input_array(2*i-1)) + 1
             nelem(input_array(2*i))   = nelem(input_array(2*i)) + 1
          ENDDO

          iopt=ppm_param_alloc_fit
          DO i = 1,nvertices
             ldc=nelem(i)
             CALL ppm_alloc(edges_per_node(i)%adj_edge,ldc,iopt,info)
             or_fail_alloc("Could not allocate edges_per_node%adj_edge")
          ENDDO

          end_subroutine()
          END SUBROUTINE allocate_processor_lists

          !---------------------------------------------------------------------
          !  For each processor, an adjacency list is formed, f.e. IF there
          !  exists an edge 1-5, adj. list of proc. 1 will contain 5 and
          !  adj. list of proc. 5 will contain 1.
          !---------------------------------------------------------------------
          SUBROUTINE fill_processor_lists(input_array,info)

          IMPLICIT NONE

          INTEGER, DIMENSION(:), INTENT(IN   ) :: input_array
          INTEGER,               INTENT(  OUT) :: info

          INTEGER :: i
          INTEGER :: pos

          start_subroutine("fill_processor_lists")

          ALLOCATE(offset(1:nvertices),STAT=info)
          or_fail_alloc("Could not allocate offset")

          offset = 1

          DO i = 1,nedges
             pos = input_array(2*i-1)
             edges_per_node(pos)%adj_edge(offset(pos)) = i
             offset(pos) = offset(pos) + 1

             pos = input_array(2*i)
             edges_per_node(pos)%adj_edge(offset(pos)) = i
             offset(pos) = offset(pos) + 1
          ENDDO

          end_subroutine()
          END SUBROUTINE fill_processor_lists

          !---------------------------------------------------------------------
          ! Allocates lists for adj. edges of each edge
          !---------------------------------------------------------------------
          SUBROUTINE allocate_edge_lists(input_array,info)

          IMPLICIT NONE

          INTEGER, DIMENSION(:), INTENT(IN   ) :: input_array
          INTEGER,               INTENT(  OUT) :: info

          INTEGER :: i,iopt
          INTEGER :: size1
          INTEGER :: size2
          INTEGER :: node1
          INTEGER :: node2

          start_subroutine("allocate_edge_lists")

          ALLOCATE(lists(1:nedges),STAT=info)
          or_fail_alloc("Could not allocate lists")

          iopt=ppm_param_alloc_fit
          DO i = 1, nedges
             node1 = input_array(2*i-1)
             node2 = input_array(2*i)
             size1 = SIZE(edges_per_node(node1)%adj_edge)
             size2 = SIZE(edges_per_node(node2)%adj_edge)

             ldc=size1+size2-2
             ! -2 comes from the fact that both edges will contain the
             ! edge itself. So, two edge wont take place and they occupy
             ! a place of 2 elements

             CALL ppm_alloc(lists(i)%adj_edge,ldc,iopt,info)
             or_fail_alloc("Could not allocate lists%ajd_edge")
          ENDDO

          end_subroutine()
          END SUBROUTINE allocate_edge_lists

          !---------------------------------------------------------------------
          ! Adjacent edges of edges are found, so that line graph is formed
          !---------------------------------------------------------------------
          SUBROUTINE fill_edge_lists(input_array,info)

          IMPLICIT NONE

          INTEGER, DIMENSION(:), INTENT(IN   ) :: input_array
          INTEGER,               INTENT(  OUT) :: info

          INTEGER :: i
          INTEGER :: j
          INTEGER :: idx
          INTEGER :: size1
          INTEGER :: size2
          INTEGER :: node1
          INTEGER :: node2
          INTEGER :: v1
          INTEGER :: v2

          start_subroutine("fill_edge_lists")

          DEALLOCATE(offset,STAT=info)
          or_fail_dealloc("Could not deallocate offset")

          ALLOCATE(offset(1:nedges),STAT=info)
          or_fail_alloc("Could not allocate offset")

          offset = 1 !Same offset array is used also in this subrout.

          DO i = 1, nedges
             node1 = input_array(2*i-1)
             node2 = input_array(2*i)

             size1 = SIZE(edges_per_node(node1)%adj_edge)

             DO j = 1, size1
                idx = edges_per_node(node1)%adj_edge(j)
                v1 = input_array(2*idx-1)
                v2 = input_array(2*idx)

                IF ((node1.NE.v1).OR.(node2.NE.v2)) THEN
                   lists(i)%adj_edge(offset(i)) = idx
                   offset(i) = offset(i) + 1
                ENDIF
             ENDDO

             size2 = SIZE(edges_per_node(node2)%adj_edge)
             DO j = 1, size2
                idx = edges_per_node(node2)%adj_edge(j)
                v1 = input_array(2*idx-1)
                v2 = input_array(2*idx)

                IF ((node1.NE.v1).OR.(node2.NE.v2)) THEN
                   lists(i)%adj_edge(offset(i)) = idx
                   offset(i) = offset(i) + 1
                ENDIF
             ENDDO
          ENDDO

          end_subroutine()
          END SUBROUTINE fill_edge_lists

          !---------------------------------------------------------------------
          ! Adjacency lists are assigned to vertices of the line graph
          !---------------------------------------------------------------------
          SUBROUTINE assign_edge_lists(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i

          start_subroutine("assign_edge_lists")

          ALLOCATE(node(1:nedges),STAT=info)
          or_fail_alloc("Could not allocate node")

          DO i=1,nedges
             node(i)%list => lists(i)%adj_edge
          ENDDO

          end_subroutine()
          END SUBROUTINE assign_edge_lists

          !---------------------------------------------------------------------
          ! Sets max_degree to delta + 1 as this is the number of colors
          ! that will be used
          !---------------------------------------------------------------------
          SUBROUTINE get_maximum_degree(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i

          start_subroutine("get_maximum_degree")

          max_degree = 0
          DO i = 1, nedges
             IF (SIZE(node(i)%list).GT.max_degree) THEN
                max_degree = SIZE(node(i)%list)
             ENDIF
          ENDDO

          max_degree = max_degree + 1

          end_subroutine()
          END SUBROUTINE get_maximum_degree

          !---------------------------------------------------------------------
          !  initialization of heap list
          !---------------------------------------------------------------------
          SUBROUTINE initialize_binary_heap(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("initialize_binary_heap")

          ALLOCATE(node_sat(0:max_degree,1:nedges),STAT=info)
          or_fail_alloc("Could not allocate node_sat")

          ALLOCATE(size_heap(0:max_degree),STAT=info)
          or_fail_alloc("Could not allocate size_heap")

          node_sat  = -1
          size_heap = 0

          end_subroutine()
          END SUBROUTINE initialize_binary_heap

          !---------------------------------------------------------------------
          !  initializes nodes
          !---------------------------------------------------------------------
          SUBROUTINE initialize_nodes(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i

          start_subroutine("initialize_nodes")

          ncolor = max_degree                       ! colors is max_degree

          DO i=1,nedges
             node(i)%color     = 0                  ! no color
             node(i)%iscolored = .FALSE.            ! noone is colored yet
             node(i)%degree    = SIZE(node(i)%list) ! number of adjacent elements
             node(i)%loc_heap  = 0                  ! location on heap list is 0
             node(i)%dsat      = 0                  ! dsat-value is 0 for all
          ENDDO

          ALLOCATE(used_color(0:ncolor),STAT=info)
          or_fail_alloc("Could not allocate used_color")

          used_color = .FALSE.

          end_subroutine()
          END SUBROUTINE initialize_nodes

          !---------------------------------------------------------------------
          ! Given the index of the node, this subroutine inserts the node in the
          ! corresponding row of binary heap lists
          !---------------------------------------------------------------------
          SUBROUTINE insert_node(idx,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: loc_parent
          INTEGER :: dsat_value
          INTEGER :: degree_own
          INTEGER :: degree_parent
          INTEGER :: loc_own

          start_subroutine("insert_node")

          dsat_value                                 = node(idx)%dsat
          size_heap(dsat_value)                      = size_heap(dsat_value) + 1
          node_sat(dsat_value,size_heap(dsat_value)) = idx
          node(idx)%loc_heap                         = size_heap(dsat_value)
          degree_own                                 = node(idx)%degree
          loc_own                                    = node(idx)%loc_heap

     !    after inserting, heapIFication must be done
          DO WHILE(loc_own.GT.1)
             loc_parent = loc_own/2
             degree_parent = node(node_sat(dsat_value,loc_parent))%degree

             IF (degree_own .GT. degree_parent) THEN
                CALL swap_nodes(node_sat(dsat_value,loc_parent),node_sat(dsat_value,loc_own),info)
                or_fail("swap_nodes")
             ELSE
                RETURN
             ENDIF
             loc_own = node(idx)%loc_heap
          ENDDO

          end_subroutine()
          END SUBROUTINE insert_node

          !---------------------------------------------------------------------
          !  swaps nodes given their IDs
          !---------------------------------------------------------------------
          SUBROUTINE swap_nodes(idx1,idx2,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx1
          INTEGER, INTENT(IN   ) :: idx2
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: temp
          INTEGER :: dsat_value

          start_subroutine("swap_nodes")

          dsat_value = node(idx1)%dsat
          !swap elements
          temp = node_sat(dsat_value,node(idx1)%loc_heap)
          node_sat(dsat_value,node(idx1)%loc_heap) = node_sat(dsat_value,node(idx2)%loc_heap)
          node_sat(dsat_value,node(idx2)%loc_heap) = temp

          temp = node(idx1)%loc_heap
          node(idx1)%loc_heap = node(idx2)%loc_heap
          node(idx2)%loc_heap = temp

          end_subroutine()
          END SUBROUTINE swap_nodes

          !---------------------------------------------------------------------
          !  Gets first element of heap list where dsat is the greatest,
          !  such that the node to be colored fulfills the max. degree
          !  among those with max. dsat value condition
          !---------------------------------------------------------------------
          FUNCTION next_node_to_color() RESULT(idx)

          IMPLICIT NONE

          INTEGER :: idx

          INTEGER :: i

          start_function("next_node_to_color")

          ! Look from top to bottom, and return the first node of the greatest
          ! dsat-value row in binary heap lists
          DO i = max_degree, 0, -1
             IF (size_heap(i) .GT. 0) THEN
                idx = node_sat(i, 1)
                GOTO 9999
             ENDIF
          ENDDO

          end_function()
          END FUNCTION next_node_to_color

          !---------------------------------------------------------------------
          !  Given the index number of the node, colors the node with
          !  minimum available color
          !---------------------------------------------------------------------
          SUBROUTINE color_edge(idx,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: color_idx
          INTEGER :: i

          start_subroutine("color_edge")

          ! initialize used_color array to FALSE
          DO i = 1,ncolor
             used_color(i) = .FALSE.
          ENDDO

          ! set used_color arrays elements to TRUE if that color was used
          DO i = 1,SIZE(node(idx)%list)
             color_idx = node(idx)%list(i)
             used_color(node(color_idx)%color) = .TRUE.
          ENDDO

          ! get the minimum color that was not used and color the edge with it
          DO i = 1, ncolor
             IF (.NOT.used_color(i)) THEN
                node(idx)%color = i
                node(idx)%iscolored = .TRUE.
                GOTO 9999
             ENDIF
          ENDDO

          end_subroutine()
          END SUBROUTINE color_edge

          !---------------------------------------------------------------------
          !  Removes the colored node from binary heap list and updates
          !  its neighbors, such that if the dsat value of the neighbor
          !  has changed, it is removed from the heaplist and inserted in
          !  the new list
          !---------------------------------------------------------------------
          SUBROUTINE update_neighbors(idx,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx !index
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i
          INTEGER :: neighbor

          start_subroutine("update_neighbors")

          ! delete the node as it was already colored
          CALL delete_node(idx,info)
          or_fail("delete_node")

          DO i = 1,SIZE(node(idx)%list)
             neighbor = node(idx)%list(i) ! get index of neighbor
             CALL compute_dsat(neighbor,info)  ! compute neighbors dsat-values
             or_fail("compute_dsat")
          ENDDO

          end_subroutine()
          END SUBROUTINE update_neighbors

          !---------------------------------------------------------------------
          !  Given the index of the node, removes it from the binary list
          !---------------------------------------------------------------------
          SUBROUTINE delete_node(idx,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx !index
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: dsat_value
          INTEGER :: idx_last

          start_subroutine("delete_node")

          dsat_value = node(idx)%dsat     !sat value of node

          IF (size_heap(dsat_value).GT.0) THEN
             node(idx)%degree = -1   !degree is set to -1 to minimize
             idx_last = node_sat(dsat_value, size_heap(dsat_value))
             CALL swap_nodes(idx,idx_last,info)
             or_fail("swap_nodes")

             CALL max_heapify(idx,info)   !the node is pushed back
             or_fail("max_heapify")

             size_heap(dsat_value) = size_heap(dsat_value) - 1 !size decreases
             node(idx)%degree = SIZE(node(idx)%list)
          ENDIF

          end_subroutine()
          END SUBROUTINE delete_node

          !---------------------------------------------------------------------
          !  Computes dsat-value of the node that index is provided as input
          !  and if dsat-value has changed, removes it from binary heap list
          !  then, inserts it back in the new row that it is supposed to be
          !---------------------------------------------------------------------
          SUBROUTINE compute_dsat(idx,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx !index
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i
          INTEGER :: neighbor
          INTEGER :: dsat_new

          start_subroutine("compute_dsat")

          used_color = .FALSE.

          ! if the node is already colored set its dsat-value to -1 and return
          IF (node(idx)%iscolored) THEN
             node(idx)%dsat = -1
          ! if the node was not colored, compute its dsat-value and update in
          ! binary heap list, if necessary
          ELSE
             ! set used_colors elements to TRUE if that color was used
             ! in at least one of the neighbors.
             DO i = 1,SIZE(node(idx)%list)
                neighbor = node(idx)%list(i)
                used_color(node(neighbor)%color) = .TRUE.
             ENDDO

             ! compute dsat-value of the node
             dsat_new = 0
             DO i = 1, ncolor
                IF (used_color(i)) THEN
                   dsat_new = dsat_new + 1
                ENDIF
             ENDDO

             ! if the computed dsat-value is greater than previous value
             IF (dsat_new .GT. node(idx)%dsat) THEN
                CALL delete_node(idx,info)     ! remove node from binary heap list
                or_fail("delete_node")

                node(idx)%dsat = dsat_new ! set dsat-value to computed one
                CALL insert_node(idx,info)     ! insert the node back in the list
                or_fail("insert_node")
             ENDIF
          ENDIF

          end_subroutine()
          END SUBROUTINE compute_dsat

          !---------------------------------------------------------------------
          ! Recursively heapifies the heap list such that top node is max.
          !---------------------------------------------------------------------
          RECURSIVE SUBROUTINE max_heapify(idx,info)

          IMPLICIT NONE

          INTEGER, INTENT(IN   ) :: idx
          INTEGER, INTENT(  OUT) :: info

          INTEGER :: loc_largest
          INTEGER :: dsat_value
          INTEGER :: loc

          start_subroutine("max_heapify")

          dsat_value  = node(idx)%dsat
          loc         = node(idx)%loc_heap
          loc_largest = loc

          IF (size_heap(dsat_value).GE.(2*loc)) THEN
             IF (node(node_sat(dsat_value,2*loc))%degree.GT. &
             &   node(node_sat(dsat_value,loc))%degree) THEN
                loc_largest = 2*loc
             ELSE
                loc_largest = loc
             ENDIF
          ENDIF

          IF (size_heap(dsat_value).GE.(2*loc+1)) THEN
             IF (node(node_sat(dsat_value,2*loc+1))%degree.GT. &
             &   node(node_sat(dsat_value,loc_largest))%degree) THEN
                loc_largest = 2*loc + 1
             ENDIF
          ENDIF

          IF (loc_largest.NE.loc) THEN !IF swapping necessary
             CALL swap_nodes(node_sat(dsat_value,loc_largest),node_sat(dsat_value,loc),info)
             or_fail("swap_nodes")

             !propogating max_heap property
             CALL max_heapify(node_sat(dsat_value,loc_largest),info)
             or_fail("max_heapify")
          ENDIF

          end_subroutine()
          END SUBROUTINE max_heapify

      END SUBROUTINE ppm_color_edge
