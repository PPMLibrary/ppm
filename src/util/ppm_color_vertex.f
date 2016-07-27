      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_color_vertex
      !-------------------------------------------------------------------------
      ! Copyright (c) 2016 MOSAIC Group (MPI-CBG Dresden)
      !
      !
      ! This file is part of the PPM_RC program.
      !
      ! PPM_RC is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM_RC is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM_RC. If not,
      ! see <http://www.gnu.org/licenses/>.
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_color_vertex(numV,nneigh,neigh,coloring,info)
        !!! This algorithm is to find coloring of a graph
        !!!
        !!! Input and Output arrays look like :
        !!!
        !!!  1           2         3         4     ...  For Vertices 1:numV
        !!!
        !!!  3,          2,        4,        0     ...  nneigh
        !!!  2,10,11,   1,3,   2,15,18,19,         ...   neigh
        !!!  1,          2,        1         1     ...  coloring
        !!!
        !!!  Reference:
        !!!  Thomas F. Coleman and Jorge J. More, Estimation of sparse Jacobian
        !!!  matrices and graph coloring problems. J. Numer. Anal. V20, P187-209, 1983

        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        IMPLICIT NONE

        INTEGER,                 INTENT(IN   ) :: numV
        !!! Number if vertices
        INTEGER, DIMENSION(:),   INTENT(IN   ) :: nneigh
        !!! Number of neighbors for every vertices
        INTEGER, DIMENSION(:,:), INTENT(IN   ) :: neigh
        !!! Neighbors for every vertices
        INTEGER, DIMENSION(:),   INTENT(  OUT) :: coloring
        !!! The color of the vertex i will be stored in coloring(i).
        !!! i.e., vertex i belongs to coloring coloring(i)
        INTEGER,                 INTENT(  OUT) :: info
        !!! Status
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(numV) :: mark
        INTEGER                  :: max_color
        INTEGER                  :: i,j,k

        CHARACTER(LEN=ppm_char) :: caller='ppm_color_vertex'

        CALL substart(caller,t0,info)

        IF (numV.LT.1) THEN
           fail("There should be at least one vertex!",ppm_error=ppm_error_fatal)
        ENDIF

        max_color=1
        mark=ppm_big_i
        coloring=numV

        !the vertices are ordered
        DO i=1,numV
           DO j=1,nneigh(i)
              k=coloring(neigh(j,i))
              mark(k)=i
           ENDDO !j=1,nneigh(i)

           j=1
           DO WHILE (j.LE.max_color.AND.mark(j).EQ.i)
              j=j+1
           ENDDO

           ! At this point, j is the smallest possible color
           coloring(i)=j

           ! All colors are used up. Add one more color
           IF (j.EQ.max_color) max_color=max_color+1
        ENDDO !i=1,numV

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_color_vertex